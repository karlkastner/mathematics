% Thu  6 Jan 15:52:37 CET 2022
fc = 3;
fun  = @(x,fc,p) (x./(x.^2 + fc.^2)).^(2*p);
%Ifun = @(fc,p) 1./(4*gamma(p-1)./gamma((p-1)/2).^2*fc^(p-1))
Ifun = @scale_bandpass_continuous;
tol = 1e-10;

dp = 0.1;
p = [0.5:0.25:4];
%2+dp:dp:4;
I = [];
for idx=1:length(p)
	a = tol.^(1/(2*p(idx)))
	xlim = [-((-(2*a*fc - 1)*(2*a*fc + 1))^(1/2) - 1)/(2*a)
		 ((-(2*a*fc - 1)*(2*a*fc + 1))^(1/2) + 1)/(2*a)]

	I(idx,1) = Ifun(fc,p(idx));
	I(idx,2) = quad(@(x) fun(x,fc,p(idx)),xlim(1),xlim(2))
	
end



