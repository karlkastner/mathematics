% Tue 31 Jan 11:14:24 CET 2023
fc=logspace(-1,1,10);
Sc=logspace(-1,1,10)';
[a,b]=lognpdf_mode2par(fc,Sc);
lognpdf(fc,a,b)-Sc
