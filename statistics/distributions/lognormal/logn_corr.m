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
%% function corr_eaeb = logn_corr(lr,lmu_a,lmu_b,lsd_a,lsd_b)
%%
%% correlation of two log-normal random variables, where the log of the variables
%% is correlated with correlation r
function corr_eaeb = logn_corr(lr,lmu_a,lmu_b,lsd_a,lsd_b)
	if (nargin()<2)
		lmu_a = 0;
	end
	if (nargin()<3)
		lmu_b = 0;
	end
	if (nargin()<4)
		lsd_a = 1;
	end
	if (nargin()<5)
		lsd_b = 1;
	end
	% standard deviation
	sd_ea = logn_std(lmu_a,lsd_a);
	sd_eb = logn_std(lmu_b,lsd_b);
	% covariance
	cov_eaeb = logn_cov(lr,lmu_a,lmu_b,lsd_a,lsd_b);
	% correlation
	corr_eaeb = cov_eaeb ./(sd_ea.*sd_eb);
end

