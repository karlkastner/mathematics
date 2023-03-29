% Mon 23 Jan 17:07:31 CET 2017
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
		A = vander_1d([-dt; 0; dt],2);

		% TODO extrapolate end-points linearly
		val_ = [val(1); val; val(end)];
		% this has shifted by 1 due to padding one value
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

