% Tue 13 Mar 12:20:15 CET 2018
% Karl Kastner, Berlin
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
%
%% finite difference derivative matrix in two dimensions
% function [Dx,Dy,D2x,Dxy,D2y] = derivative_matrix_2d(n,L,order,bc)
function [Dx,Dy,D2x,Dxy,D2y] = derivative_matrix_2d(n,L,order,bc)
	if (1 == length(order))
		order(2) = order(1);
	end
	% first oder derivative matrices
	Dx1  = derivative_matrix_1_1d(n(1),L(1),order,bc{1});
	Dy1  = derivative_matrix_1_1d(n(2),L(2),order,bc{2});

	% bug corrected on
	% Fri 18 May 10:04:14 CEST 2018
	Dx  = kron(speye(n(2)),Dx1);
	Dy  = kron(Dy1,speye(n(1)));

	if (abs(order)<2)
		order = 2;
	end

	% second order derivative
	D2x1 = derivative_matrix_2_1d(n(1),L(1),order,bc{1});
	D2y1 = derivative_matrix_2_1d(n(2),L(2),order,bc{2});

	% this yields still a 1-neighbourhood kernel with 9 points of which 4
	% are non-zero, as the differentiation is on orthogonal directions
	%Dxy = kron(Dy1,Dx1); % This was a bug!
	% nb : this is identical to Dy*Dx
	Dxy = Dx*Dy;

	% order corrected (swapped) on Wed 29 Aug 12:19:24 CEST 2018
	D2x = kron(speye(n(2)),D2x1);
	D2y = kron(D2y1,speye(n(1)));
end

