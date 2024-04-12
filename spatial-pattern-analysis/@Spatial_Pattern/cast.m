% Mon  4 Dec 17:45:13 CET 2023
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
function obj = cast(obj,type)
	field_C = {'x','y','r'};
	for idx=1:length(field_C)
		obj.(field_C{idx}) = type(obj.(field_C{idx}));
	end
	obj.f = cast_deep(obj.f,@single);
	obj.w = cast_deep(obj.w,@single);
	obj.S = cast_deep(obj.S,@single);
end

