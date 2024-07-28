x = linspace(0,100,1e4)';
fx = fourier_axis(x);
fc_ = [];
n=10;
Sc=logspace(-1,1,10);
 fc=logspace(-1,1,10);
 for idx=1:n;
 for jdx=1:10;
 [a,b]=lognpdf_mode2par(fc(idx),Sc(jdx));
 Sc_exact(idx,jdx) = Sc(jdx);
 fc_exact(idx,jdx) = fc(idx);
 S = lognpdf(fx,a,b);
 [Sc__(idx,jdx),mdx] = max(S);
 fc__(idx,jdx) = abs(fx(mdx));
 [fc_(idx,jdx), Sc_(idx,jdx)] = lognpdf_mode(a,b);
 end;
 end;
 Sc_exact
 Sc_
 Sc__
 fc_exact
 fc_
 fc__
