% Wed 29 Mar 16:49:51 CEST 2023
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
function lpar = geometric_ar1_2d_grid_cell_averaged_moment2par(mu,sd,theta,dx)
	% TODO for use analytic solution for dx=0 as initial value
	lpar = [1,1];
	lpar = lsqnonlin(resfun,lpar);

	function res = objective(lpar)
		mu_ = exp(lpar(1) + 0.5*lpar(2));
		sd_ = geometric_ar1_2d_grid_cell_averaged_std(lpar(1),lpar(2),theta,dx);
		res = [mu_-mu,sd_-sd]; 
		
	end
end

