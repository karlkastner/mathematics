% Thu 26 Mar 12:46:07 +08 2020
n=1e4; x=lognrnd(2,3,n,1); J.func=@geomean; J.estimate(x); J.serr, std(x)/sqrt(n-1), geoserr(x)
