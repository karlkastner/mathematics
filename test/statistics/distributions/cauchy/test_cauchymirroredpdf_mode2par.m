% 2024-06-29 16:35:09.175788999 +0200
fc=1;
Sc=1;
fx = linspace(0,3);
[a,b]=cauchymirroredpdf_mode2par(fc,Sc)
[fc,Sc]=cauchymirroredpdf_mode(a,b)
S =cauchymirroredpdf(fx,a,b);
plot(fx,S)

