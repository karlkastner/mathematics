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
%
%% determine covariance between grid cell averged values of a stationary
%% stochastic process on an equispaced grid
%%
%% cov(e_ij,e_kl) = E((e_ij - mu)*(e_kl - mu))
%%                = E(e_ij e_kl) - mu
%%                = E(int int e dx1 dy1 int int e dx2 dy2) - mu
%%                = E(int int int e e dx1 dy1 dx2 dy2) - mu
%%                = E(int int int s^2 rho(x2-x1,y2-y1) dx1 dy1 dx2 dy2) - mu
%%
%% f_ij = int_(x_i - dx/2)^(x_i+dx/2) int_(y_i-dy/2)^(y_j+dy/2) f(x,y) dx dy
%%
%% integrals approximated by Gauss method
%%
function cov_ = cov_grid_cel_averages(rfun,mu,s,x,y,order)
	dx = x(2)-x(1);
	dy = y(2)-y(1);
	xx = flat(repmat(cvec(x),1,length(y)));
	yy = flat(repmat(rvec(y),length(x),1));
	% weights and evaluation points of gauss integral
	[w,p] = int_1d_gauss(order);
	% integrate along x for first cell
	x1 = 0   + p*(0.5*dx*[-1;1])
	% integrate along x for second cell
	x2 = xx + (p*(0.5*dx*[-1;1]))';
	% integrate along y for first cell
	y1 = 0   + p*(0.5*dy*[-1;1])
	% integrate along y for second cell
	y2 = yy + (p*(0.5*dy*[-1;1]))';
	cov_ = zeros(numel(xx),1);
	% 2D gauss integration
	for x1id=1:length(w)
	 for x2id=1:length(w)
          for y1id=1:length(w)
           for y2id=1:length(w)
		% note it is not necessary to multiply or divide by dx and dy,
		% as they the multiplication for integration cancels with the division for averaging
		cov_ = cov_ + w(x1id)*w(x2id)*w(y1id)*w(y2id)*rfun(x1(x1id)-x2(:,x2id),y1(y1id)-y2(:,y2id));
		%cov_ = cov_ + w(x1id)*w(x2id)*w(y1id)*w(y2id)*rfun(x1(x1id),x2(:,x2id),y1(y1id),y2(:,y2id));
	   end % for y2id
	  end % for yi1d
	 end % for x2id
	end % for x1id
	cov_ = sd*sd*reshape(cov_,length(x),length(y)) - mu;
end % cov_

