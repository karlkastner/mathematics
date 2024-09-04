% 2023-02-13 13:57:28.894060160 +0100

fm = 1.5;
s = 0.9;
f = linspace(0,3);
[fc,Sc] = normalmirroredpdf_mode(fm,s);
S = normalmirroredpdf(f,fm,s);
[fm_,s_] = normalmirroredpdf_mode2par(fc,Sc)
S_ = normalmirroredpdf(f,fm_,s_);

clf;
plot(f,S);

hold on
plot(fc,Sc,'*')

plot(f,S_)

