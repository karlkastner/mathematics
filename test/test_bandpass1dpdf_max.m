
% low-pass
f0 = [1,2,3];
p  = [1,2,3];
%fun = @(fr,f0,p) 1./(1 + (fr/f0).^2).^(2*p);
L= 100;
n = L^2;
f = linspace(0,L,n)';


'int S'
for idx=1:length(f0)
 for jdx=1:length(p)
	%Sc(idx,jdx) = integral(@(f) lowpass1dpdf(f,f0(idx),p(jdx),0),0,inf);
	% max occurs always at 0
	S = bandpass1dpdf(f,f0(idx),p(jdx),true);
	[Sc(idx,jdx),mdx] = max(S);
	fc(idx,jdx) = f(mdx);

	%[fc_(idx,jdx),Sc_(idx,jdx)] = bandpass1dpdf_mode(f0(idx),p(jdx));
	[Sc_(idx,jdx)] = bandpass1dpdf_max(f0(idx),p(jdx));
 end
end
Sc
Sc_
max(abs(Sc-Sc_),[],'all')

% gamma(p+2)*gamma(p-1)/(gamma(2*p))/(2*p+2)

