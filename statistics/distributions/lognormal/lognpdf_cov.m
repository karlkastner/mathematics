% Wed 29 Mar 11:12:38 CEST 2023
% Karl Kastner, Berlin
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%% covariance of two log-normally distributed random variables,
%
%% cov(ea,eb) = cov(exp(mua + sa*za),exp(mub + sb*zb))
%% where za, zb are standard normal distributed and correlated
function cov_eaeb = lognpdf_cov(mu_za, mu_zb, s_za, s_zb, r)
	if (nargin()<2)
		mu_za = 0;
	end
	if (nargin()<3)
		mu_zb = 0;
	end
	if (nargin()<4)
		s_za = 1;
	end
	if (nargin()<5)
		s_zb = 1;
	end
	% cov(ea,eb) = E( (ea-Eea)(ea - Eeb)^2 )
	%            = E(ea eb) - Eea Eeb
	% Eea        = E[exp(za)] = exp(mua + 1/2 sa^2)
	% E(ea*eb)   = exp(mua + mub) E exp(sa*za + sb*zb)
	% split into two uncorrelated
	%	     = exp(mua + mub) * E exp(sa za + sb zb)
	%	     = exp(mua + mub) * E exp(sa za + sb (r za + sqrt(1-r^2) zc))
	% by independence
	%	     = exp(mua + mub) * E exp((sa + r*sb) za) E exp( sb sqrt(1-r^2) zc)
	% by expectation of the log-normal distribution
	%	     = exp(mua + mub) * exp((sa + r*sb)^2/2) exp(sb^2 (1-r^2)/2)
	%	     = exp(mua + mub) * exp((sa^2 + 2 r sa sb + r^2*sb^2)/2 + sb^2 (1-r^2)/2)
	%	     = exp(mua + mub) * exp((sa^2 + 2 r sa sb + sb^2)/2)
		%	     = exp(1 + r

	% this neglects the contribution by mu_za, scaled in later
	Eea      = lognpdf_mean(mu_za,s_za);
	Eeb      = lognpdf_mean(mu_zb,s_zb);
	Eeaeb    = exp(0.5*(s_za.*s_za + 2*r*s_za.*s_zb + s_zb.*s_zb));
	cov_eaeb = exp(mu_za + mu_zb).*Eeaeb - Eea*Eeb;
end

