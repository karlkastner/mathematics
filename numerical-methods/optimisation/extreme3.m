% Mon 23 Jan 17:07:31 CET 2017
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
%% extract maxima by quadratic approximation from sampled function val(t)
%% intended to be called after [mval, mid] = max(val) for refinement of
%% locatian and maximum
%%
%% input
%% t    : sampling time (uniformly spaced)
%% v    : values at sampling times
%% ouput:
%% tdx  : index where extremum should be computed
%% t0   : location of the extremum
%% val0 : value of extremum
%%
% TODO take t0 and dt as arguments
% TODO no need for two step approach with max,
%      just compute extremum for each double-interval
% TODO automatic switch between uniformly and not-uniformly sampled data
% function [val0, t0, ddv_dt2] = extreme3(t,val,tdx)
function [val0, t0, ddv_dt2] = extreme3(t,val,tdx)
		val = cvec(val);
		t   = cvec(t);
		nt  = length(t);
		dt  = t(2)-t(1);

		% quadratic vandermonde matrix
		% TODO allow for changing time step
		t_ = [t(1); t; t(end)];
		% shift local to reduce round of error
		t3 = [t_(tdx)-t_(tdx+1); zeros(1,length(tdx)); t_(tdx+2)-t_(tdx+1)];
		% A = vander_1d([-dt; 0; dt],2);
		A = vander_1d(t3,2);

		% TODO extrapolate end-points linearly
		val_ = [val(1); val; val(end)];
		% the shift of 1 compensates for the padding of one value
		val3 = [val_(tdx); val_(tdx+1); val_(tdx+2)];
		%val3 = [val(max(1,tdx-1)); val(tdx); val(min(nt,tdx+1))];

		% TODO matrix mul is associative, just precompute pd*inv(A)
		% fit a quadratic polymomial
		c0   = A \ double(val3);
		% derive
		% c1 = polyder_man(fipud(c0));
		c1  = [c0(2,:); 2*c0(3,:)];

		% find roots (zeros)
		dt0 = -c1(1,:)./c1(2,:);

		% limit
		fdx = abs(dt0) > dt;
		dt0(fdx) = 0;

		%% v'(dt0) = 0 and v''(dt0) determines type of extremum
		ddv_dt2 = 2*c0(2,:);

		% maximum value
		val0 = c0(1,:) + c0(2,:).*dt0 + c0(3,:).*dt0.^2;
		t0   = t(tdx) + dt0;
end

