% Mon 13 Feb 13:42:27 CET 2023
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
function [fm,sf] = normpdf_wrapped_mode2par(fc,Sc)
	% type requirement for lsqnonlin
	fc = double(fc);
	Sc = double(Sc);

	par0 = [fc,0.1./Sc];
	% [par,resnorm,residual,exitflag] 
	par = lsqnonlin(@resfun,par0);
	fm = par(1);
	sf = par(2);
	
function res = resfun(par)
	[fc_,Sc_] = normpdf_wrapped_mode(par(1),par(2));
	res = [fc_-fc,Sc_-Sc];
end

end

