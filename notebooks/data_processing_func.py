from datetime import date, datetime, time
from typing import Optional
import pandas as pd
import numpy as np
if __name__ == '__main__':
    raise RuntimeError("You should import this file instead of run it directly")

calendar_map = {
    "2025-2026": {
        "sem1": {"start": "2025-09-01", "end": "2025-11-29"},
        "sem2": {"start": "2026-01-05", "end": "2026-04-18"},
        "summer": {"start": "2026-05-11", "end": "2026-06-27"}
    },
    "2026-2027": {
        "sem1": {"start": "2026-09-07", "end": "2026-12-19"},
        "sem2": {"start": "2027-01-11", "end": "2027-04-18"},
        "summer": {"start": "2027-05-17", "end": "2027-07-03"}
    }
}

def _opening_fraction(index, overrides):
    values = []
    for t in index:
        day = t.strftime("%Y-%m-%d")
        hours = overrides.get(day, ["08:20", "22:00"] if t.weekday() < 5
                              else ["08:20", "19:00"] if t.weekday() == 5
                              else ["11:00", "19:00"])
        if hours is None:
            values.append(0.0)
            continue
        start, end = time.fromisoformat(hours[0]), time.fromisoformat(hours[1]) # why do not have error
        if start >= end:
            raise ValueError("Opening overrides must start before closing (no overnight hours)")
        opening = t.normalize() + pd.Timedelta(hours=start.hour, minutes=start.minute)
        closing = t.normalize() + pd.Timedelta(hours=end.hour, minutes=end.minute)
        overlap = (min(t + pd.Timedelta(hours=1), closing) - max(t, opening)).total_seconds()
        values.append(max(0.0, overlap) / 3600)
    return np.asarray(values)


def _acf(values, lag):
    z = values - values.mean()
    denominator = float(np.dot(z, z))
    if denominator <= 1e-12:
        return 1.0 if lag == 0 else 0.0
    return 1.0 if lag == 0 else float(np.dot(z[:-lag], z[lag:]) / denominator)


def _features(y, overrides, calendar = calendar_map):

    x = pd.DataFrame(index=y.index)
    for lag in (1, 24, 168):
        x[f"lag_{lag}"] = y.shift(lag)
    past = y.shift(1)
    x["trend_6"] = past.rolling(6, min_periods=6).mean()
    hour = y.index.hour
    weekly_hour = y.index.dayofweek * 24 + hour
    for k in range(1, 6):
        x[f"daily_sin_{k}"] = np.sin(2*np.pi*k*hour/24)
        x[f"daily_cos_{k}"] = np.cos(2*np.pi*k*hour/24)
        x[f"weekly_sin_{k}"] = np.sin(2*np.pi*k*weekly_hour/168)
        x[f"weekly_cos_{k}"] = np.cos(2*np.pi*k*weekly_hour/168)
    flags = np.zeros(len(y), dtype=np.int8)
    dates = y.index.strftime("%Y-%m-%d")
    for academic_year in calendar.values():
        for key, flag in (("sem1", 1), ("sem2", 2), ("summer", 3)):
            if key in academic_year:
                span = academic_year[key]
                flags[(dates >= span["start"]) & (dates <= span["end"])] = flag
    x["semester_flag"] = flags
    x["weekday"] = y.index.dayofweek
    x["opening_fraction"] = _opening_fraction(y.index, overrides)
    x["available_last_168"] = past.rolling(168, min_periods=1).count()
    for lag in range(11):
        x[f"acf_lag_{lag}"] = past.rolling(168, min_periods=168).apply(
            lambda values, n=lag: _acf(values, n), raw=True)
    return x.astype(float)


def data_reshape(data: dict|pd.DataFrame,
                 timedelta_hours: Optional[str] = "1h"):
#=============================================================================
  # cols = ["recorded_at", "occupant_change", "occupant_count"]
  time_format = "%Y-%m-%d_%H:%M:%S"
  if type(data) == "dict":
    result = pd.DataFrame(data,
                  index=data["recorded_at"],
                  columns=["occupant_change", "occupant_count"]
                )
  else:
    result = data.copy()
    result = result.set_index("recorded_at")
#=============================================================================
  result.index = pd.to_datetime(result.index, format=time_format)

  result['Entry'] = np.where(result['occupant_change'] == 1, result['occupant_change'], 0)
  result['Exit'] = np.where(result['occupant_change'] == -1, result['occupant_change'], 0)

  result.drop(columns=['occupant_change'], inplace=True)

  result_without_count = result.drop(columns=['occupant_count'])
  result_count = result.loc[:, 'occupant_count']

  #set(result_without_count.index.duplicated())
  # Remark: cause someone entry in sametime is possible
  result_without_count = result_without_count.groupby('recorded_at', as_index=True).sum()
  result_count = result_count.groupby('recorded_at', as_index=True).mean()

  result_resample_woc = result_without_count.resample(timedelta_hours).sum()
  result_resample_count = result_count.resample(timedelta_hours).mean()

  result_resample_woc = result_resample_woc.between_time('08:20', '22:00')
  result_resample_count = result_resample_count.between_time('08:20', '22:00')

#=============================================================================

  return result_resample_woc, result_resample_count
