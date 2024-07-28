fc=1;Sc=1;fx = linspace(0,3); [a,b]=cauchywrappedpdf_mode2par(fc,Sc), [fc,Sc]=cauchywrappedpdf_mode(a,b), S =cauchywrappedpdf(fx,a,b); plot(fx,S)

