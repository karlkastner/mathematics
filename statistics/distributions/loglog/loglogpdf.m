% 2019-02-19 20:02:26.459675962 +0100
%
% PDF of the log-log distribution
%
function f = loglogpdf(x,mu,s)
	%  x > 1
	f = 1./(sqrt(2*pi*s)*x.*log(x)).*exp(-log(log(x))-mu).^2./(2*s));
end

