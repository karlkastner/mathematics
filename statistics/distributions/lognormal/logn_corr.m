% Wed 29 Mar 11:20:01 CEST 2023
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
%% correlation of two log-normal random variables, where the log of the variables
%% is correlated with r
function corr_eaeb = corr_logn(r,mu_za,mu_zb,s_za,s_zb)
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
	% standard deviation
	s_ea = sqrt(exp(s_za.*s_za) - 1).*exp(mu_za + 0.5*s_za.*s_za);
	s_eb = sqrt(exp(s_zb.*s_zb) - 1).*exp(mu_zb + 0.5*s_zb.*s_zb);
	% covariance
	cov_eaeb = logn_cov(r,mu_za,mu_zb,s_za,s_zb);
	% correlation
	corr_eaeb = cov_eaeb ./(s_ea.*s_eb);
end

