% Tue  7 Nov 10:43:18 CET 2023
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
% circular boundary conditions
% TODO other boundary conditions
function A = upsampling_matrix(n)
	%A = 2*downsampling_matrix(2*n)';
	%w = [1/6,2/3,1/6];
	w = [1/2,1,1/2];
	A = spdiags(ones(n,1)*w,-1:1,n,n);
	A = A(1:2:end,:);
	A(1,end) = 0.5;
	A = A';
%	A = [3/4;1/4];
%	B = kron(eye(n),1-A);
%	A = kron(eye(n),A),
%	A = A+circshift(B,[0,+1]);
%	A = circshift(A,[-1,-1]);
end

