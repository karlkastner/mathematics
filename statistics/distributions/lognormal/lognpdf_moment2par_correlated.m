% Wed 29 Mar 16:23:12 CEST 2023
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
%% mean, variance and correlation of the log of two correlated log-normal random varibales
%%
%% let
%%	z_i = exp(lz_i),	i = 1,2
%% with
%%	E(z_i) = mui
%%	Var(z_i) = sd_i^2
%%	corr(z_1,z_2) = r
%%
%% then this function determines the log-moments:
%%
%% 	E(lz_i) = lmu_i
%%	Var(lz_i) = lsd_i^2
%%	corr(lz_1,lz)2) = lr
%%
%% which are identical with the parameters of the log-normal functions in matlab
%%
function [lmu1,lmu2,lsd1,lsd2,lr] = lognpdf_moment2par_correlated(mu1,mu2,sd1,sd2,r)
	[lmu1,lsd1]     = lognpdf_moment2par(mu1,sd1);
	[lmu2,lsd2]     = lognpdf_moment2par(mu2,sd2);
	% lr = log(1 + r.*(exp(lsd.^2) - 1))./lsd.^2;
	lr = log(r.*sqrt((exp(lsd1.^2) - 1).*(exp(lsd2.^2) - 1)) + 1)./(lsd1.*lsd2);
end

