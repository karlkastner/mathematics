% Fri 18 Oct 19:52:45 CEST 2024
%
%% Kalman-filter
classdef Kalman < handle
	properties
		x
		% transition matrix
		A
		H
		R
		Q
		P
	end % properties
		
	methods
		function obj = Kalman()
		end

		function x_k = update(obj,z)
	
			% last state
			x_km1 = obj.x;
	
			% predict the new state after transition
			% x = A*x + w
			x_k = obj.A*x_km1;
		
			% z : measurement
			% x : states (aka parameter of the RD-model)
			% z = H*x
	
			% project the covariance
			P = obj.A*obj.P*obj.A' + obj.Q;
				
			% correct estimates
			% Kalman gain
			K = P*obj.H' / (obj.H*P*obj.H' + obj.R);
	
			% update estimate
			x_k = x_k + K*(z - obj.H*x_k);
	
			% update covariance
			P = (eye(length(obj.x)) - K*obj.H)*P;
		
			% write back
			obj.x = x_k;
		end % function update
	end % methods
end % classdef Kalman

