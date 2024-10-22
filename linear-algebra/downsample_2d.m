% 2024-01-04 19:52:41.137057243 +0100
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
%% note: hw is twice as fast as fw
%% x = rand(2^10); tic; for idx=1:100; x_=x; for jdx=1:10; x_=downsample_2d(x_,'fw'); end;end; toc
function x = downsample_2d(x,mode)
	switch (mode)
	case {'fw'}
		x = downsample1(x);
		x = downsample2(x);
	case {'hw'}
		x = (  0.5*x(1:2:end-1,1:2:end-1) ... % centre
		    + 0.125*x(1:2:end-1,2:2:end) ... % right
		    + 0.125*[x(1:2:end-1,end),     x(1:2:end-1,2:2:end-2)] ... % lef
		    + 0.125*x(2:2:end,1:2:end)  ... % down
		    + 0.125*[x(end,1:2:end-1);     x(2:2:end-2,1:2:end-1)] ... % up
		    );
	end
end

