% 2015-12-17 12:27:51.575538568 +0100
%% Kalman filter

function X = kalman(z,u)
	if (nargin() < 2)
		u = zeros(size(z));
	end
	
	% prediction matrix (constant mean)
	F = 1;
	% transformation to observed to internal state
	H = 1;
	% transformation of control input to state
	B = 1;
	% observation covariance
	R = 0.1;
	% initial covariance matrix
	P = 1;
	
	% process noise
	% how to estimate?
	% large Q means low confidence in model prediction
	% Q = FR_vF', R_v = (v-\bar v)(v-\bar v)' (sensor noise)
	% F = dx/dv
%	Q = ?
	Q = 0.1;

	% initial state
	x = z(1);

	X = zeros(size(z));
	X(1) = z(1);
	for idx=2:length(z)
		% predict state
		x = F*x + B*u(:,idx); % + v (unknown noise)
		X(idx) = x;
		% predict covariance
		P = F*P*F' + Q;
		% measurement residual
		y = z(:,idx) - H*x;
		% residual covariance
		S = H*P*H' + R
		% optimal Kalman gain
		K = P*H' * inv(S);
		% correct estimate
		x = x + K*y;
		% update of the covariance matrix
		I = 1;
		P = (I-K*H)*P;
	end
end % kalman

