% Wed 29 Mar 21:26:05 CEST 2023
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%
%% determine covariance between grid cell averged values of a stationary
%% stochastic process on an equispaced grid
%%
%% cov(e_ij,e_kl) = E[ (e_ij - mu)*(e_kl - mu) ]
%%                = 1/dx^2 E[ int (e(x1) - mu) dx1 int (e(x2)-mu) dx2 ]
%%                = 1/dx^2 E[ int int (e(x1) - mu) (e(x2) - mu) dx1 dx2 ) ]
%%                = 1/dx^2 int int cov(x2-x1) dx1 dx2
%%
%% f_ij = int_(x_i - dx/2)^(x_i+dx/2) f(x) dx
%%
%% integrals approximated by Gauss' method
%%
function cov_ = cov_cel_averages_1d(cfun,x,dx,m,efficient)
	if (nargin()<4)
		m = 2;
	end
	if (nargin()<5)
		efficient = 1;
	end
	siz = size(x);
	x = flat(x);
	% weights and evaluation points of gauss integral

if (efficient)
	dxs  = dx/m;
	cov_ = m*cfun(x);
	for idx=1:m-1
		% note that +/- dx is only identical when x = 0
		cov_ = cov_ + (m-idx)*cfun(x+dxs*idx);
		cov_ = cov_ + (m-idx)*cfun(x-dxs*idx);
	end
	cov_ = cov_ / (m*m);
else
	[w,p] = int_1d_gauss(order);
	%[w,p] = int_1d_equal(order);

	% integrate along x for first cell
	x1 = 0   + p*(0.5*dx*[-1;1]);
	% integrate along x for second cell
	x2 = x + (p*(0.5*dx*[-1;1]))';

	cov_ = zeros(prod(siz),1);
	% 1D gauss integration
	for x1id=1:length(w)
	 for x2id=1:length(w)
		% note it is not necessary to multiply or divide by dx,
		% as they the multiplication for integration cancels with the division for averaging
		cov_ = cov_ + w(x1id)*w(x2id)*cfun(x1(x1id)-x2(:,x2id));
	 end % for x2id
	end % for x1id
	cov_ = reshape(cov_,siz(1),siz(2));
end
end % cov_

