% Tue 28 Feb 14:45:20 CET 2023
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
% ouput:
%	cov_    : covariance bettween grid cell averages, for cells averaged
%		  at x,y and 0,0
%
% input
%	cfun	: covariance function of the stationary process
%	x,y	: grid cell centre coordinates, to be correlated with cell around 0,0
%	m,n	: subgrid cell number
%
%% determine covariance between grid cell averged values of a stationary
%% stochastic process on an equispaced grid
%%
%% cov(e_ij,e_kl) = E[ (e_ij - mu)*(e_kl - mu) ]
%%		  = E[ (e_00 - mu)*(e_(k-i,l-j) - mu]
%%                = 1/dx^2 1/dy^2 E[ int int (e(x1,y1) - mu) dx1 dy1 int int (e(x2,y2)-mu) dx2 dy2 ]
%%                = 1/dx^2 1/dy^2 E[ int int int int (e(x1,y1) - mu) (e(x2,y2) - mu) dx1 dy1 dx2 dy2) ]
%%                = 1/dx^2 1/dy^2   int int int int cov(x2-x1,y2-y1) dx1 dy1 dx2 dy2
%%
%% f_ij = int_(x_i - dx/2)^(x_i+dx/2) int_(y_i-dy/2)^(y_j+dy/2) f(x,y) dx dy
%%
%% integrals approximated by equal spaced mid-point intervals,
%% this allows to reduce the double-integral along each dimension into a
%% single integral and hence to reduce the computational effort from m^4 to m^2
%%
% note it is not necessary to multiply or divide by dx and dy,
% as they the multiplication for integration cancels with the division for averaging
function cov_ = cov_cell_averages_2d(cfun,x,y,dx,dy,m,n,efficient)
	if (nargin()<6)
		m = 1;
	end
	if (nargin()<7)
		n = m;
	end
	if (nargin()<8)
		efficient = 1;
	end

	if (efficient)
		% efficient integration, double-integrals along each dimension reduced
		% to single integrals by exploiting stationarity
		% sub-interval-length
		dxs = dx/m;
		dys = dy/n;

		% no shift
		cov_ = m*n*cfun(x,y);
		for idx=1:m-1
			% shift only one coordinate
			cov_ = cov_ + (m-idx)*n*cfun(x-dxs*idx,y);
			cov_ = cov_ + (m-idx)*n*cfun(x+dxs*idx,y);
			cov_ = cov_ + m*(n-idx)*cfun(x        ,y-dys*idx);
			cov_ = cov_ + m*(n-idx)*cfun(x        ,y+dys*idx);
			for jdx=1:m-1
				% shift both coordinates
				cov_ = cov_ + (m-idx)*(n-jdx)*cfun(x-dxs*idx,y-dys*jdx);
				cov_ = cov_ + (m-idx)*(n-jdx)*cfun(x-dxs*idx,y+dys*jdx);
				cov_ = cov_ + (m-idx)*(n-jdx)*cfun(x+dxs*idx,y-dys*jdx);
				cov_ = cov_ + (m-idx)*(n-jdx)*cfun(x+dxs*idx,y+dys*jdx);
			end % for jdx
		end
		cov_ = cov_/(m*m*n*n);
	else
		% explicitly solve the double integral (for testing)
		siz = size(x);
	
		x = flat(x);
		y = flat(y);
		% weights and evaluation points of gauss integral
		%[w,p] = int_1d_gauss(order);
		[w,p] = int_1d_equal(m);
	
		% integrate along x for first cell
		x1 = 0   + p*(0.5*dx*[-1;1]);
		% integrate along x for second cell
		x2 = x + (p*(0.5*dx*[-1;1]))';
		% integrate along y for first cell
		y1 = 0   + p*(0.5*dy*[-1;1]);
		% integrate along y for second cell
		y2 = y + (p*(0.5*dy*[-1;1]))';
		cov_ = zeros(prod(siz),1);
	
		% 2D gauss integration
		for x1id=1:length(w)
		 for x2id=1:length(w)
	          for y1id=1:length(w)
	           for y2id=1:length(w)
			cov_ = cov_ + w(x1id)*w(x2id)*w(y1id)*w(y2id)*cfun(x1(x1id)-x2(:,x2id),y1(y1id)-y2(:,y2id));
			%cov_ = cov_ + w(x1id)*w(x2id)*w(y1id)*w(y2id)*rfun(x1(x1id),x2(:,x2id),y1(y1id),y2(:,y2id));
		   end % for y2id
		  end % for yi1d
		 end % for x2id
		end % for x1id
		cov_ = reshape(cov_,siz(1),siz(2));
	end
end % cov_cell_averages_2d

