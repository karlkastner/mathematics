fm = 1.1;
s = 0.9;
f = linspace(0,3);
S = normpdf_wrapped(f,fm,s);
clf;
plot(f,S);
[fc,Sc] = normpdf_wrapped_mode(fm,s);
hold on
plot(fc,Sc,'*')

[fm_,s_] = normpdf_wrapped_mode2par(fc,Sc)
S = normpdf_wrapped(f,fm_,s_);
plot(f,S)

