import pandas as pd
import numpy as np

from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.utils.validation import check_X_y, check_array, check_is_fitted
import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, RegressorMixin
from statsforecast import StatsForecast

class Weekly_mean(BaseEstimator, RegressorMixin):
    def __init__(self):
        self.weekday_means_ = None

    def fit(self, y: pd.Series):
        if not isinstance(y.index, pd.DatetimeIndex):
            raise ValueError("y must have a DatetimeIndex")
        g = y.groupby(y.index.weekday).mean()
        arr = np.full(7, float(y.mean()), dtype=float)   # fill missing weekdays with overall mean
        for k, v in g.items():
            arr[int(k)] = float(v)
        self.weekday_means_ = arr
        return self

    def predict(self, dates):
        check_is_fitted(self, "weekday_means_")
        idx = dates if isinstance(dates, pd.DatetimeIndex) else pd.to_datetime(pd.Index(dates))
        return self.weekday_means_[idx.weekday.astype(int)]

# statsforecast MODELS
"""
Single variable models 
AutoETS
AutoTheta
AutoTBATS
AutoCES


DynamicOptimizedTheta

baseline
HistoricAverage
Naive
RandomWalkWithDrift
SeasonalNaive
WindowAverage
SeasonalWindowAverage

target +  variable models 
AutoMFLES	
"""

from statsforecast.models import (
    # Auto method
    AutoETS,
    AutoTheta,
    AutoTBATS,
    AutoCES,
    AutoARIMA,
    # Theta
    DynamicOptimizedTheta,
    # Baseline
    HistoricAverage,
    Naive,
    RandomWalkWithDrift,
    SeasonalNaive,
    WindowAverage,
    SeasonalWindowAverage
)

season_length_weekly = 13
window_size_3_periods = 3
models = [
    # Automatic exponential smoothing (ETS)
    AutoETS(),  # statsforecast will try select ETS automatically

    # Theta-family models
    AutoTheta(),  # automatic Theta selection
    DynamicOptimizedTheta(),  # optimized/dynamic theta variant

    # TBATS-like automatic model (can be slow / heavy)
    AutoTBATS(
        use_boxcox =True,
        use_trend=True,
        use_damped_trend=True,
        use_arma_errors=False,  # toggle depending on your needs
        season_length=[24, season_length_weekly],  # candidate seasonal periods
    ),

    # Complex exponential smoothing / damped exponential smoothing family
    AutoCES(),  # automatic CES configuration

    AutoARIMA(),

    # --- Baselines ---
    HistoricAverage(),  # mean over history
    Naive(),  # last-observed value repeated
    RandomWalkWithDrift(),  # random walk with linear drift
    SeasonalNaive(season_length=season_length_weekly),  # repeat value from same hour-of-week
    WindowAverage(window_size= window_size_3_periods),  # average over last 7 days (adjust window)
    SeasonalWindowAverage(window_size= window_size_3_periods, season_length=season_length_weekly),  # seasonal windowed average
]

class SFEstimator(BaseEstimator, RegressorMixin):
    """
    Lightweight adapter: single-series sklearn-style estimator that uses StatsForecast internally.
    Efficient: builds a StatsForecast instance in fit() and reuses it for predict().
    """
    def __init__(self, models, freq='H', fallback_model=None, n_jobs=1):
        """
        models: list of statsforecast model instances (e.g., [SeasonalNaive(season_length=168)])
        freq: pandas frequency string for the series
        fallback_model: a model instance to use if needed (passed to StatsForecast)
        n_jobs: passed to StatsForecast
        """
        #const
        self.models = models
        self.freq = freq
        self.fallback_model = fallback_model
        self.n_jobs = n_jobs
        #var
        self._sf = None
        self.train_series_ = None
        self.df_ = None

    def fit(self, y: pd.Series):
        if not isinstance(y.index, pd.DatetimeIndex):
            raise ValueError("y must be a pandas Series with a DatetimeIndex")
        
        self._sf = StatsForecast(models=self.models, freq=self.freq, fallback_model=self.fallback_model, n_jobs=self.n_jobs)
        self.train_series_ = y.sort_index()

        # build panel DataFrame expected by StatsForecast
        self.df_  = pd.DataFrame({"unique_id": "series_0", "ds": self.train_series_.index, "y": self.train_series_.values})
        print(self.df_.tail(5))
        self._sf.fit(self.df_ )
        return self

    def predict(self, h_or_dates):
        """
        fcst: ['unique_id','ds', '<modelname>']
        Warning:
        StatsForecast may accept ds directly depending on version
        """
        check_is_fitted(self, "_sf")
        # determine horizon and target dates
        if isinstance(h_or_dates, (int, np.integer)):
            h = int(h_or_dates)
            fcst = self._sf.forecast(h=h, df = self.df_)
            # fcst columns: ['unique_id','ds', '<modelname>']
            return fcst
        else:
            # accept list-like of timestamps or pd.DatetimeIndex
            ds = h_or_dates if isinstance(h_or_dates, pd.DatetimeIndex) else pd.to_datetime(pd.Index(h_or_dates))
            h = len(ds)
            fcst = self._sf.forecast(h=h, df = self.df_)
            return fcst
