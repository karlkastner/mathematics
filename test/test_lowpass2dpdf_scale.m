
% low-pass
f0 = [1,2,3];
p  = [1,2,3];
%fun = @(fr,f0,p) 1./(1 + (fr/f0).^2).^(2*p);

for idx=1:length(f0)
 for jdx=1:length(p)
	Sc(idx,jdx) = integral(@(f) 2*pi*f.*lowpass2dpdf(0,f,f0(idx),p(jdx),true),0,inf);
 end
end
Sc
max(abs(Sc-1),[],'all')

