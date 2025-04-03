% 2025-02-25 15:15:20.645614021
% Karl Kästner, Berlin
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.
%

% covariance of exp(z) is lower than covariance of z
% covariance of exp(z) decreases with sigma 

n=1e6;
 lr = 0.9
 lmu1=1;
 lmu2=1/2;
 ls1=1/3;
 ls2=1/5;
 x = randn(n,1);
 y = lr*x + sqrt(1-lr^2)*randn(n,1);
 x = lmu1 +ls1*x;
 y=lmu2+ls2*y;
 x = exp(x);
 y = exp(y);
'mu'
 mean(x)
 exp(lmu1+1/2*ls1^2)
 mean(y)
 exp(lmu2+1/2*ls2^2)
 'var'
 var(x)
% var_ = exp(2*mu1+2*s1^2) - exp(mu1+1/2*s1^2)^2
 var1 = exp(2*lmu1+ls1^2)*(exp(ls1^2) - 1)
 var(y)
 var2 = exp(2*lmu2+ls2^2)*(exp(ls2^2) - 1)
 'cov'
 cov(x,y)
% cov_ = exp(z1)*exp(z2)	- exp(mu1+1/2*s1^2)*exp(mu2+1/2*s2^2)
% cov_ = exp(z1+z2)	- exp(mu1+1/2*s1^2)*exp(mu2+1/2*s2^2)
% cov_ = exp(mu1 + mu2 + e1 + e2)	- exp(mu1+mu2+1/2*s1^2+1/2*s2^2)
 cov_ = exp(mu1 + mu2 + lr*s1*s2 + 1/2*(s1^2+s2^2)) - exp(mu1+mu2+1/2*s1^2+1/2*s2^2)
 cov_ = exp(mu1 + mu2 + 1/2*(s1^2+s2^2))*(exp(lr*s1*s2) - 1)
% cov_ = exp(2*mu+(1+r)*s^2) - exp(mu+1/2*s^2)^2
% cov_ = exp(2*mu+s^2)*(exp(r*s^2) - 1)
 'corr'
 corr(x,y)
 corr12 = cov_/sqrt(var1*var2)
% cov_ = exp(mu1 + mu2 + 1/2*(s1^2+s2^2))*(exp(r*s1*s2) - 1) ...
%	/ sqrt(exp(2*mu1+s1^2)*(exp(s1^2) - 1)*exp(2*mu2+s2^2)*(exp(s2^2) - 1))
 corr12_ = (exp(lr*s1*s2) - 1) ...
	/ sqrt((exp(s1^2) - 1)*(exp(s2^2) - 1))
% cov_ = (exp(r*s^2) - 1)/(exp(s^2) - 1)
 cfun = @(r,s) (exp(r.*s.^2) - 1)./(exp(s.^2) - 1)
% corr_=  (exp(2*mu+(1+r)*s^2) - exp(mu+1/2*s^2)^2)
x = linspace(-50,50,1e4);
 r = exp(-abs(x))';
 c=[r,cfun(r,[0,0.5,1])];
subplot(2,2,1)
 plot(x,c);
subplot(2,2,2)
 S=abs(fft(c)).^2;
plot(S./S(1,:))
xlim([0,80]);

syms lm ls x t a; r = exp(-x/t); c=lognpdf_corr(lm,lm,ls,ls,r); c=simplify(c); s=solve(c==str2sym('exp(-1)')), s=simplify(s); sfun = matlabFunction(s); sfun(10,100)
 
