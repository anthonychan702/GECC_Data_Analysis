###-------------------------------------------------------------###
### STAT 3005 
### Chapter 4
###-------------------------------------------------------------###

###-------------------------------------------------------------###
### This file contains general functions used in Chapter 4
###-------------------------------------------------------------###
### (Part I) Sample correlation
###   - general:         cor0 <--------------- Recommended 
###   - pearson:         cor.p0, cor.p1
###   - spearman:        cor.s0, cor.s1, cor.s2, cor.s3
###   - kendall:         cor.k0, cor.k1, cor.k2, cor.k3
###   - bergsma-dassios: cov.bd0, cov.bd1, cov.bd2 
###   - chatterjee:      cor.c0, cor.c1
###
### (Part II) Indepedence test 
###   - general:         cor.test0 <--------------- Recommended 
###   - pearson:         cor.test.p0, cor.test.p1
###   - spearman:        cor.test.s0, cor.test.s1 
###   - kendall:         cor.test.k0, cor.test.k1 
###   - bergsma-dassios: cov.test.bd0 
###   - chatterjee:      cor.test.c0, cor.test.c1
###
###-------------------------------------------------------------###
 







###------------------------------------------------------------------
### (Part I) Sample correlation
###------------------------------------------------------------------
###------------------------------------
### General 
###------------------------------------
cor0 = function(x,y,method=c("pearson","spearman","kendall","bergsma-dassios","chatterjee")){
	method = match.arg(method)
	switch(method,
		"pearson"=cor.p0(x,y),
		"spearman"=cor.s0(x,y),
		"kendall"=cor.k0(x,y),
		"bergsma-dassios"=cov.bd0(x,y,uStatistic=TRUE),
		"chatterjee"=cor.c0(x,y))
}

###------------------------------------
### Pearson's correlation 
###------------------------------------
cor.p0 = function(x,y){  
	cor(x,y)
}
cor.p1 = function(x,y){   
	(mean(x*y)-mean(x)*mean(y))/sqrt((mean(x^2)-mean(x)^2)*(mean(y^2)-mean(y)^2))
}

###------------------------------------
### Spearman's correlation 
###------------------------------------
cor.s0 = function(x,y){  
	cor(x,y,method="spearman")
}
cor.s1 = function(x,y){   
	cor(rank(x),rank(y))
}
cor.s2 = function(x,y){    
	n = length(x)
	1-6*sum((rank(x)-rank(y))^2)/n/(n^2-1)
}
cor.s3 = function(x,y){     
	n = length(x)
	cor(1:n,rank(y[order(x)]))
}

###------------------------------------
### Kendall's correlation 
###------------------------------------
cor.k0 = function(x,y){
	cor(x,y,method="kendall")
}
cor.k1 = function(x,y){
	n = length(x)
	s = 0
	for(i in 2:n){
		for(j in 1:(i-1)){
			s = s+sign((x[i]-x[j])*(y[i]-y[j]))
		}
	}
	s/choose(n,2)
}
cor.k2 = function(x,y){
	n = length(x)
	s = 0
	for(i in 2:n){
		s = s+sum(sign((x[i]-x[1:(i-1)])*(y[i]-y[1:(i-1)])))
	}
	s/choose(n,2)
}
cor.k3 = function(x,y){
	pcaPP::cor.fk(x,y)  
}


###------------------------------------
### Bergsma--Dassios's covariance 
###------------------------------------
aFUN = function(z1,z2,z3,z4,p=1){
	temp = abs(z1-z2)^p+abs(z3-z4)^p-abs(z1-z3)^p-abs(z2-z4)^p
	sign(ifelse(abs(temp)<1e-10,0,temp))
}
cov.bd1 = function(x,y,uStatistic=TRUE){
	n = length(x)
	c = 0
	if(uStatistic){
		for(i1 in 1:n){
			for(i2 in (1:n)[-i1]){
				for(i3 in (1:n)[-c(i1,i2)]){
					c = c+sum(aFUN(x[i1],x[i2],x[i3],x[-c(i1,i2,i3)])*aFUN(y[i1],y[i2],y[i3],y[-c(i1,i2,i3)]))
				}
			}
		}
		d = ifelse(c<0, -1, 1)
		out = d*exp(log(d*c)-log(n)-log(n-1)-log(n-2)-log(n-3)) 
	}else{
		for(i1 in 1:n){
			for(i2 in 1:n){
				for(i3 in 1:n){
					c = c+sum(aFUN(x[i1],x[i2],x[i3],x)*aFUN(y[i1],y[i2],y[i3],y))
				}
			}
		}
		d = ifelse(c<0, -1, 1)
		out = d*exp(log(d*c)-4*log(n))	
	}
	out
}
cov.bd2 = function(x,y,J=100000){
	# subsample version 
	n = length(x)
	s = 0
	i = array(sample(1:n,J*4,replace=TRUE),dim=c(J,4))
	for(j in 1:J){
		s = s+sum(   aFUN(x[i[j,1]],x[i[j,2]],x[i[j,3]],x[i[j,4]])
				   * aFUN(y[i[j,1]],y[i[j,2]],y[i[j,3]],y[i[j,4]]) )
	}
	s/J
}
#-------------------------------------------------------------------------------------
# BD covariance (fast version)
# - an O(n^2) algorithm that allows ties; see Heller and Heller (2016) arXiv:1605.08732
# - modified from the cpp code retrived from the "TauStar" R pakcage 
#-------------------------------------------------------------------------------------
vecToRanks = function(x){
    n = length(x)
    xSortedInds = order(x)-1 
    ranks = repeats = rep(0,n)
    k = 1 
    last = x[xSortedInds[1]+1]
    for(i in 0:(n-1)){
        current = x[xSortedInds[i+1]+1] 
        if(current != last){
            k = k+1
            last = current 
        }
        repeats[i+1] = k
    }
    for(i in 0:(n-1)) {
        ranks[xSortedInds[i+1]+1] = repeats[i+1] 
    }
    ranks
}
ranksToLeqMat = function(xRanks,yRanks){
	n = length(xRanks)
	xMax = max(xRanks)
	yMax = max(yRanks)
    leqMat = array(0, dim=c(xMax+1,yMax+1))
    for(i in 0:(n-1)) {
        leqMat[xRanks[i+1]+1, yRanks[i+1]+1] = leqMat[xRanks[i+1]+1, yRanks[i+1]+1] + 1
    }
    for(i in 1:xMax){
        for(j in 1:yMax){
            leqMat[i+1,j+1] = leqMat[i+1,j] + leqMat[i,j+1] - leqMat[i, j] + leqMat[i+1,j+1]
        }
    }
    leqMat
}
leqMatToUniqueCountMat = function(leqMat){
	xMax = nrow(leqMat)-1
	yMax = ncol(leqMat)-1
    uCountMat = array(0, dim=c(xMax+1,yMax+1))
    for(i in 1:(xMax)) {
        for(j in 1:(yMax)){
            numEqYAndLessOrEqX = leqMat[i+1, j+1] - leqMat[i+1, j]
            uCountMat[i+1,j+1] = uCountMat[i+1, j] + (numEqYAndLessOrEqX * (numEqYAndLessOrEqX - 1))/2
        }
    }   
    uCountMat   
}
indexUvec = function(x, inds) {
	n = length(inds)
	newVec = rep(0,n)
	for(i in 0:(n-1)) {
		newVec[i+1] = x[inds[i+1]+1]
  	}
	newVec
}
cov.bd0 = function(x,y,uStatistic=TRUE){
	n = length(x)
	xRanks = vecToRanks(x)
	yRanks = vecToRanks(y)
	xMax = max(xRanks)
	yMax = max(yRanks)
	leqMat = ranksToLeqMat(xRanks, yRanks)
	xOrder = order(xRanks)-1
	xRanks = indexUvec(xRanks,xOrder)
	yRanks = indexUvec(yRanks,xOrder)
	uCountMat = leqMatToUniqueCountMat(leqMat)
	numCon = numDis = 0
	if(uStatistic){
		for(i in 0:(n-2)){
		    for(j in (i+1):(n-1)){
		        xRankMin = xRanks[i+1]
		        yRankMin = min(yRanks[i+1], yRanks[j+1])
		        yRankMax = max(yRanks[i+1], yRanks[j+1])
		        bot = leqMat[xRankMin, yRankMin]
		        mid = ifelse(yRankMin == yRankMax, 0, leqMat[xRankMin, yRankMax] - leqMat[xRankMin, yRankMin+1])
		        top = leqMat[xRankMin, yMax+1] - leqMat[xRankMin, yRankMax+1]
		        eqMin = leqMat[xRankMin, yRankMin+1] - leqMat[xRankMin, yRankMin]
		        eqMax = leqMat[xRankMin, yRankMax+1] - leqMat[xRankMin, yRankMax]
		        numCon = numCon + choose(top,2) + choose(bot,2)
		        if(yRankMin != yRankMax){
		            numDis = numDis + top*(mid+eqMin+bot) + bot*(mid+eqMax) + eqMin*(mid+eqMax) + eqMax*mid + choose(mid,2)
		            numDis = numDis - uCountMat[xRankMin, yRankMax] + uCountMat[xRankMin, yRankMin+1]
		        }
		    }
		}
		c = 16*numCon - 8*numDis
		d = ifelse(c<0, -1, 1)
		out = d*exp(log(d*c)-log(n)-log(n-1)-log(n-2)-log(n-3)) 
	}else{
		for(i in 0:(n-1)){
			top = leqMat[xRanks[i+1], yMax+1] - leqMat[xRanks[i+1], yRanks[i+1]+1]
			bot = leqMat[xRanks[i+1], yRanks[i+1]]
			numCon = numCon + (choose(top,2)+choose(bot,2))/2 + (top + bot)/4 
			if(i<(n-1)){
				for(j in (i+1):(n-1)){ 
					xRankMin = xRanks[i+1]
					yRankMin = min(yRanks[i+1], yRanks[j+1])
					yRankMax = max(yRanks[i+1], yRanks[j+1])
					bot = leqMat[xRankMin, yRankMin]
					mid = ifelse(yRankMin==yRankMax, 0, leqMat[xRankMin,yRankMax]-leqMat[xRankMin, yRankMin+1])
					top = leqMat[xRankMin, yMax+1] - leqMat[xRankMin, yRankMax+1]
					eqMin = leqMat[xRankMin, yRankMin+1] - leqMat[xRankMin, yRankMin]
					eqMax = leqMat[xRankMin, yRankMax+1] - leqMat[xRankMin, yRankMax]
					numCon = numCon + choose(top,2) + choose(bot,2) + (top+bot)/2 
					if(yRankMin!=yRankMax){
						numDis = numDis + top*(mid+eqMin+bot) + bot*(mid+eqMax) + eqMin*(mid+eqMax) + eqMax*mid + choose(mid,2) 
						numDis = numDis - uCountMat[xRankMin, yRankMax] + uCountMat[xRankMin, yRankMin+1]
					}
				}
			}
		}
		c = 16*numCon - 8*numDis
		d = ifelse(c<0, -1, 1)
		out = d*exp(log(d*c)-4*log(n)) 
	}
	out
}


###------------------------------------
### Chatterjee's correlation  
###------------------------------------
cor.c1 = function(x,y){
	r = rank(y[order(x)])
	n = length(x)
	1- 3*sum(abs(r[-1]-r[-n]))/(n^2-1)  
}
cor.c0 = function(x,y){
	XICOR::xicor(x,y)
}













###------------------------------------------------------------------
### (Part II) Independence test
###------------------------------------------------------------------
###------------------------------------
### General 
###------------------------------------
cor.test0 = function(x,y,method=c("pearson","spearman","kendall","bergsma-dassios","chatterjee")){
	method = match.arg(method)
	switch(method,
		"pearson"=cor.test.p0(x,y),
		"spearman"=cor.test.s0(x,y),
		"kendall"=cor.test.k0(x,y),
		"bergsma-dassios"=cov.test.bd0(x,y,nRep=1000),
		"chatterjee"=cor.test.c0(x,y))
}

###------------------------------------
### Pearson's correlation 
###------------------------------------
cor.test.p0 = function(x,y){  
	cor.test(x,y)$p.value
}
cor.test.p1 = function(x,y){   
	n = length(x)
	cor = cor(x,y)
	1-pchisq(n*cor^2,1)
}

###------------------------------------
### Spearman's correlation test
###------------------------------------
cor.test.s0 = function(x,y){
	cor.test(x,y,method="s")$p.value
}
cor.test.s1 = function(x,y){ 
	n = length(x)
	cor = cor(rank(x),rank(y))  
	t = cor/sqrt(var)  
	p = 1-pchisq(t^2,1)  
	p
}


###------------------------------------
### Kendall's correlation test
###------------------------------------
cor.test.k0 = function(x,y){
	cor.test(x,y,method="k")$p.value
}
cor.test.k1 = function(x,y){
	n = length(x) 
	s=0
	for(i in 2:n){
		s = s+sum(sign((x[i]-x[1:(i-1)])*(y[i]-y[1:(i-1)])))
	}
	cor = s/choose(n,2) 
	var = 2*(2*n+5)/(9*n*(n-1)) 
	p = 1-pchisq(t^2,1)  
	p
}



###------------------------------------
### Bergsma--Dassios's covariance test
###------------------------------------
cov.test.bd1 = function(x,y,uStatistic=TRUE,nRep=1000){
    n = length(x)
    out = rep(NA,nRep)
    cat(1:10,"\n")
    n0 = 50 
    for(iRep in 1:nRep){
        x0 = runif(n0)
        y0 = runif(n0)       
        out[iRep] = n0*cov.bd0(x0,y0,uStatistic=uStatistic)
        if(iRep%%floor(nRep/10)==0) cat("> ")
    }
    cat("\n")
    t = n*cov.bd0(x,y,uStatistic=uStatistic)
    p = mean(out>=t)
    p
}

#-------------------------------------------------------------------------------------
# BD covariance (U-statistic) test (fast version)
# - based on the asymptoptic distrbution when x & y ~ discrete/continous distrbution (Nandy, Weihs and Drton, 2016)
# - modified from the cpp code retrived from the "TauStar" R pakcage 
#-------------------------------------------------------------------------------------
cov.test.bd0 = function(x,y,nRep=1000){
	n = length(x)
	p1 = as.numeric(table(x))/n
	p2 = as.numeric(table(y))/n
	xIsCts = (length(p1)==n)
	yIsCts = (length(p2)==n)
	if(xIsCts&yIsCts) out = rbdcov_continuous(nRep=nRep, n0=50)
	if(xIsCts&(!yIsCts)) out = rbdcov_mixed(nRep=nRep, p1=p2, n0=50)
	if((!xIsCts)&yIsCts) out = rbdcov_mixed(nRep=nRep, p1=p1, n0=50)
	if((!xIsCts)&(!yIsCts)) out = rbdcov_discrete(nRep=nRep, p1=p1, p2=p2, n0=50)	
	cat("\nComputing the observed sample BD-covariance\n")
	t = n*cov.bd0(x,y)
	p = mean(out>=t)
	p 
}

eigenForDiscreteProbs = function(p){
	n = length(p)
	cdf = rep(NA,n)
	q1 = rep(NA,n)
	cdf[1] = p[1]
	for(i in 2:n){
		cdf[i] = cdf[i-1] + p[i]
	}
	q1[1] = p[1]*(1-cdf[1])
	for(i in 2:n){
		q1[i] = q1[i-1] + p[i] * (1-cdf[i])
	}
	symMat = array(0, dim=c(n,n))
	for(i in 1:n) {
		for(j in i:n) {
			if(i==1){
				symMat[i,j] = (1-cdf[j])^2
			}else{
				symMat[i,j] = (1-cdf[j])^2 + cdf[i-1]^2
			}
			if(i!=j) {
				symMat[i,j] = symMat[i,j] - cdf[i]*(1-cdf[i]) + (q1[j-1]-q1[i])
				symMat[j,i] = symMat[i,j]
				symMat[j,i] = symMat[j,i]*sqrt(p[i]*p[j])
			}
			symMat[i,j] = symMat[i,j]*sqrt(p[i]*p[j])
		}
		if(i%%floor(n/10)==0) cat(".")
	}
	eigen(symMat)$values
}
rbdcov_continuous = function(nRep=1000, n0=50){
	# X~continous and Y~continous
	# The vector out stores simulated values generated from the asymptoptic distribition of n*BDcorUStat
	cat("\nSimulating sample from asymptotic null distribution:\n")
	out = rep(NA,nRep)
	cat(1:10,"\n")
	for(iRep in 1:nRep){
		out[iRep] = 36*sum((rchisq(n0^2,1)-1)/rep(1:n0,each=n0)/rep(1:n0,n0))/pi^4
		if(iRep%%floor(nRep/10)==0) cat("> ")
	}
	cat("\n")
	out
}
rbdcov_discrete = function(nRep=1000, p1, p2, n0=50){
	# X~discrete and Y~discrete
	# The vector out stores simulated values generated from the asymptoptic distribition of n*BDcorUStat
	cat("Computing eiegn-values associated with the x-variable")
	eigen1 = eigenForDiscreteProbs(p1)
	cat("\nComputing eiegn-values associated with the y-variable")
    eigen2 = eigenForDiscreteProbs(p2)
    cat("\n")
    n1 = min(n0,length(eigen1))
    n2 = min(n0,length(eigen2))
    eigen1 = eigen1[1:n1]
    eigen2 = eigen1[1:n2]
    cat("\nSimulating sample from asymptotic null distribution:\n")
    out = rep(NA,nRep)
	cat(1:10,"\n")
	for(iRep in 1:nRep){
		out[iRep] = 4*sum((rchisq(n1*n2,1)-1)*rep(eigen1,each=n2)*rep(eigen2,n1)) 
		if(iRep%%floor(nRep/10)==0) cat("> ")
	}
	cat("\n")
	out
}
rbdcov_mixed = function(nRep=1000, p1, n0=50){
	# X~discrete and Y~continous
	# The vector out stores simulated values generated from the asymptoptic distribition of n*BDcorUStat
	cat("Computing eiegn-values associated with the discrete variable")
	eigen1 = eigenForDiscreteProbs(p1)
	n1 = min(n0,length(eigen1))
    eigen1 = eigen1[1:n1]
    cat("\nSimulating sample from asymptotic null distribution:\n")
    out = rep(NA,nRep)
	cat(1:10,"\n")
	for(iRep in 1:nRep){
		out[iRep] = 12*sum((rchisq(n1*n0,1)-1)*rep(eigen1,each=n0)/rep(1:n0,n1))/pi^2
		if(iRep%%floor(nRep/10)==0) cat("> ")
	}
	cat("\n")
	out
}

###------------------------------------
### Chatterjee's correlation test
###------------------------------------
cor.test.c1 = function(x,y){
	r = rank(y[order(x)])
	n = length(x)
	cor = 1- 3*sum(abs(r[-1]-r[-n]))/(n^2-1) 
	var = 2/5/n 
	t = cor/sqrt(var) 	
	p = 1-pnorm(t) 
	p
}
cor.test.c0 = function(x,y){
	XICOR::xicor(x,y,pvalue=TRUE)$pval
}







