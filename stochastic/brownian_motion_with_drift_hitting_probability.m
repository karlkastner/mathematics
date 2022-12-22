% Abundo 2012
% karatzas shreve 2019 p192 eq 5.12
a=1; b = 1; f = @(x) abs(a)./x.^(3/2).*normpdf((a-b*x)/sqrt(x)); quad(f,sqrt(eps),100)
