% Tue  4 Feb 16:33:21 CET 2025
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
function sk = binormpdf_skewness(p,mu1,mu2,sd1,sd2)
%E[(x-mu).^3] = E[x^3] - 3 mu E[x^2] - 3 E[x mu^2] + E[mu^3]
%		E( p x1^3 + (1-p) x2^3) - 3 mu*sd^2 - 3 mu^3 + 1 mu^3
%		E x1^2 = e^2 + 2 e mu + mu^2 = s^2 + mu^2
%		E x1^3 = e^3 + 3 e^2 mu + 3 mu^2 e + mu^3 = 3 s^2 mu + mu^2
%		= 0 - 3 (p s1^2 + (1-p)s2^2)*mu - 3 mu^3 + 1 mu^3
%		= -2 mu3 - 3 (p s1^2 + (1-p)s2^2)*mu
%	( p1 x1 + (1-p) x2)
	mu = binormpdf_mean(p,mu1,mu2,sd1,sd2);
	sd = binormpdf_std(p,mu1,mu2,sd1,sd2);

	sk = (  (      p *mu1*(3*sd1^2+mu1^2) ...
	          + (1-p)*mu2*(3*sd2^2+mu2^2)  ) ...
	      - 3*mu*(p*(sd1^2 + mu1^2) + (1-p)*(sd2^2 + mu2^2)) ...
	      + 2*mu^3 ...
	     )/sd^3;
%sk=	(p*(mu1 - mu2)*(p - 1)*( 2*mu1*mu2 + 2*mu1^2*p + 2*mu2^2*p - mu1^2 - mu2^2 - 3*sd1^2 + 3*sd2^2 - 4*mu1*mu2*p))/(p*sd1^2 - sd2^2*(p - 1) - p*(mu1 - mu2)^2*(p - 1))^(3/2)
%sk=	(p*(mu1 - mu2)*(p - 1)*( -(mu1 - mu2)^2 + 2*p*(mu1 - mu2)^2  - 3*sd1^2 + 3*sd2^2  ))/(p*sd1^2 - sd2^2*(p - 1) - p*(mu1 - mu2)^2*(p - 1))^(3/2)
%sk=	(p*(mu1 - mu2)*(p - 1)*( (2*p-1)*(mu1 - mu2)^2 - 3*sd1^2 + 3*sd2^2  ))/(p*sd1^2 - sd2^2*(p - 1) - p*(mu1 - mu2)^2*(p - 1))^(3/2)
%	sk = (-2*mu^3 - 3*(p*sd1^2 + (1-p)*sd2^2)*mu)/sd^3;
%	sk = mu*(-2*mu^2 - 3*(p*sd1^2 + (1-p)*sd2^2))/sd^3;
end

