% Tue 31 Oct 09:44:02 CET 2023
% Karl Kästner, Berlin
%
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
% wrapper for the zeros-function which supports half-percision
function x = zeros_man(varargin)
	if (ischar(varargin{end}))
		class_ = varargin{end};
		varargin = varargin(1:end-1);
	else
		class_ = 'double';
	end
	switch (class_)
	case {'half'}
		x = repmat(half(0),varargin{:});
	otherwise
		x  = zeros(varargin{:},class_);
	end
end

