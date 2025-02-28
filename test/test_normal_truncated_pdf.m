% TODO for some reason the expressions for mean and sk by Horrace are not working as interpreted, use penders expressions instead
https://real-statistics.com/normal-distribution/truncated-normal-distribution/

a=linspace(-2,2,10);
x = randn(1e6,1);

for idx=1:length(a)
fdx=x>a(idx);

mu(idx,1) = mean(x(fdx))
sd(idx,1) = std(x(fdx))
sk(idx,1) = skewness(x(fdx))
%Z = normcdf(a);
% s2  = 1 - a*normpdf(a)/Z - (normpdf(a)/Z).^2;
%sqrt(s2),
mu(idx,2) = normal_truncated_pdf_mean(0,1,a(idx))
sd(idx,2) = normal_truncated_pdf_std(0,1,a(idx))
sk(idx,2) = normal_truncated_pdf_skewness(0,1,-a(idx))
%, normal_truncated_pdf_std(0,1,-a)
end
