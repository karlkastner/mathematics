% Wed 25 Oct 11:04:13 CEST 2023
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
% impulse response of the discrete phase noise integration process
% (damped oscillator with phase noise)
% 	(z + a1 dz/dx + a2 d^2z/dx^2) == e
function [ir,A] = damped_oscillator_2d_discrete_ir(n,L,a1,a2,a3)
	Ix  = speye(n(1));
	Iy  = speye(n(2));
	D1x = derivative_matrix_1_1d(n(1),L(1),2,'circular','circular');	
	D1y = derivative_matrix_1_1d(n(2),L(2),1,'circular','circular');	
	D2x = derivative_matrix_2_1d(n(1),L(1),2,'circular','circular');
	D2y = derivative_matrix_2_1d(n(2),L(2),2,'circular','circular');
	D1x = kron(D1x,Iy);
	D2x = kron(D2x,Iy);
	D1y = kron(Ix,D1y);
	D2y = kron(Ix,D2y);
	% TODO compute D4 directly
	I   = speye(prod(n));
	% TODO inclusion of D2y into the bracket should yield more realistic results and make the
	%      pattern and pdf look more realistic, however, 
	%      this likely requires an explicit computation of the mixed derivatives,
	%      as the current implementation with just squaring A'A leads to instabilities
	%A = (I + a1*D1x + a2*D2x + a3*D2y);
	A = (I + a1*D1x + a2*D2x)*(I - a3*D2y);
	z = zeros(prod(n),1);
	z(1,1) = prod(L);
	if (0)
		tic
		%ir = bicgstabl(A,z,[],sum(n));
		%ir = cgs(A,z,[],2*sum(n));
		ir = gmres(A,z,[],[],sum(n));
		toc
	else
		tic
		ir = A \ z;
		toc
	end
	ir = reshape(ir,n);
end



