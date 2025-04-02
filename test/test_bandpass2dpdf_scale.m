
% low-pass
f0 = [1,2,3];
p  = [1,2,3];
%fun = @(fr,f0,p) 1./(1 + (fr/f0).^2).^(2*p);

normalize = true;
%normalize = false;
for idx=1:length(f0)
 for jdx=1:length(p)
	Sc(idx,jdx)  = integral(@(f) 2*pi*f.*bandpass2dpdf(0,f,f0(idx),p(jdx),normalize),0,inf);
%	Sc_(idx,jdx) = bandpass2dpdf_max(f0(idx),p(jdx));
 end
end
Sc
Sc_
max(abs(Sc-1),[],'all')

