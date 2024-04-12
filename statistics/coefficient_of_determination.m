% Sun 30 Sep 11:26:44 CEST 2018
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
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function r2 = coefficient_of_determination(y,yp,np,mode)
	if (nargin() < 4 || isempty(mode))
		mode = 'pearson';
	end
	if (nargin() < 3 || isempty(np))
		np = 1;
	end
	switch (lower(mode))
		case {'pearson'}
			rho = corr(y,yp,'type','Pearson');
		case {'hodges-lehmann'}
			rho = hodges_lehmann_correlation(y,yp);
		case {'kendall'}
			rho = kendall_to_pearson(corr(y,yp,'type','Kendall'));
		case {'spearman'}
			rho = spearman_to_pearson(corr(y,yp,'type','Spearman'));
		otherwise
			error('here');
	end % switch
	
	% number of samples
	% TODO allow for effective sample size
	ns = length(y);
	% correct for sample size
	r2 = (ns-1)/(ns-np)*rho^2;
end % coefficient_of_determination

