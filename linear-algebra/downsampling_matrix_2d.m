% 2024-01-04 20:22:25.391693409 +0100
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
%% function [R1,R2t,R1_,R2_] = downsampling_matrix_2d(n)
%%
%% downsampling matrices and restriction (coarsening) matrices for multigrid
%%
%% for vectorized 2d-data x = X2d(:)
%% 
%% R12*x = (R1+R2)*x = x_
%% where size(x) = [n(1)*n(2),1] and size(x_) = [n(1)*n(2)/4,1]
%% 
%% for operators:
%% R1*L*R2 = L_, where size(L) = (n(1)*n(2))^2 and size(L_) = 1/4*n(1)*n(2)
%%
%% made 'hw' : half-weighting
%%	downsamples with a 5-point stencil:
%%		[0,1,0; 1,4,1; 0,1,0]/8
%%	preserves sparsity of operators, i.e. when the original operatpr L
%%	has a 5-point kernel then the downsampled L_ has a 5-point kernel
%%	- note that the upsampling (prolongation) operator is not 4*the transpose
%%	  of the restriction operator, instead, the full-weighting (bilinear)
%%	  operator should be chosen for upsampling,
%%	  as hw-sampling does _not_ use information of x-cross points,
%%	  and only half the information from +-cross points
%%	  the resulting weights for using the transpose as upsampling
%%	  are 1 for retained points, 1/2 for + interpolated points and
%%	  0 for x-interpolated points
%%
%% mode 'fw' : full-weighting (bilinear interpolation)
%%	downsamples with a 9-point stencil:
%%		[1,2,1]'*[1,2,1]/16 = [1,2,1; 2 4 2; 1 2 1]/16
%% 	the downsampled operator has a 9-point kernel even when the original
%%	operator has a 5-point kernel
%%
function [R1_,R2t_,R1,R2t] = downsampling_matrix_2d(n,mode)
	R1  = downsampling_matrix_1d(n(1),'mg');
	R2t = downsampling_matrix_1d(n(2),'mg')';
	switch (mode)
	case {'fw'}
		R1_  = kron(R1,R2t');
		R2t_ = 4*R1_';
	case {'hw'}
		I  = speye(n);
		% injection operator
		%I1   = I(:,1:2:end)*sqrt(2);
 		%I2   = I(1:2:end,:)*sqrt(2);
		I1   = I(:,1:2:end);
 		I2   = I(1:2:end,:);
		R1_  = kron(R1,I2);
		R2t_ = 2*kron(I1,R2t);
		%R1_ = (R1__ + R2t__')/2;
		%R2t_ = (R1_ - 0*kron(I1,I2')'/2)';
		%R2t_ = kron(R1,R2t')';
	end
	%if (nargout()>2)
		%I1 = speye(n(1))
		%I2 = speye(n(2))
		%R1_ = kron(R1,I2);
		%R2_ = R1_';
	%	R12 = kron(R1,R2t');
		%R2_ = kron(I1,R2t);
	%end
end

