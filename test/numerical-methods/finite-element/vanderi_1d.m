% Fri  6 Jan 15:46:15 CET 2017
% Karl Kästner, Berlin
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
% vander monde matrix of the dth-derivative
%function V = vanderi_1d(x,order,d)
% TODO change definition to make d last argument
function V = vanderi_1d(x,d,order)
	for idx=1:order+1
		V(:,idx) = factorial(idx-1)./factorial(idx+d-1)*x.^(idx-1+d);
	end

%	if (0 == d)
%		% no derivative
%		V = vander_1d(x,order);
%	else
%		V      = zeros(size(x,1),order+1,class(x));
%		if (d<=order)
%		V(:,d+1) = factorial(d);
%	%	for idx=d+1:order+1
	%		V(:,idx) = V(:,idx-1).*x;
	%	end
	%	for idx=d:order+1
	%		V(:,idx) = factorial(idx)/factorial(idx-d)*V(:,idx);
	%	end
%		for idx=d+2:order+1
%			V(:,idx) = (idx-1)/(idx-d-1)*V(:,idx-1).*(x);
%		end
%		end
%	end
end % vanderi_1d

