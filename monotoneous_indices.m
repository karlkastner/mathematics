% Tue 31 Jan 15:27:52 CET 2023
% Karl Kästner, Berlin
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
function fdx = monotoneous_indices(x,mode)
	x= cvec(x);
	switch (mode)
	case {'descending','d',-1}
		xmin = cummin(x);
		fdx  = [true;x(2:end)<xmin(1:end-1)];
	case {'ascending','a',1}
		xmax = cummax(x);
		fdx  = [true;x(2:end)>xmax(1:end-1)];
	otherwise
		error('unknown mode')
		disp('mode');
	end
end

