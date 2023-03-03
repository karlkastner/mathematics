% Fri 24 Feb 13:31:13 CET 2023

T = 2;
nt = 1000;
t = (0:nt)'/nt*T;

S0 = 3;
r  = 0.1;
sigma = 2;

n = 100000;

S_sd = gbb_std(t,S0,t(end),r,sigma);

S    = gbb_simulate(t,r,sigma,S0,n);

%sqrt(mean(var(S(end,:))))

S_sd_ = geostd(S,[],2);
S_mu = geomean(S,2);

plot(t,[S_sd,S_sd_]);
%,S_sd./S_sd_,sqrt(S_sd)])
legend('analytic','numeric')
'entire series'
gsd = gbb_geostd_entire_series(t,S0,T,r,sigma) 
gsd(2) = sqrt(geomean(S_sd.^2))
gsd(3) = sqrt(geomean(S_sd_.^2))
gsd/gsd(1)

if (1)
syms t S0 T r sigma
gsd = gbb_std(t,S0,T,r,sigma)
%gsd=exp((sigma*t^(1/2)*(T - t)^(1/2))/T^(1/2))
gsd = simplify(gsd,'ignoreanalyticconstraints',true)
end


