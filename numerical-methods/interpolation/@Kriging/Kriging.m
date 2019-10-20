% Mon Jun 16 07:19:29 WIB 2014
% Karl Kastner, Berlin
%% class for Kriging interpolation
% TODO make Order and Model sub-classes of Kriging
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
classdef Kriging < Interpolator

	properties
		% parameter of set of the semivariance model
		param;
		% standard deviation of the estimated semivariance model parameter
		sparam;
		cparam;
		lincparam;
		varhist
		% threshold for invoking iterative solver
		MAX_SIZE_DIRECT = 1000;
		% tolerance for iterative solver
		% the uncertainty in the data is high, no need to solve with high accuracy
		TOL = 1e-3;
		% cuf off threshhold, measured points x_in where the estimated variance
		% with respect to the interpolated point x_out reaches 100*tresh per-cent
		% of the global variance will be cut off
		thresh = 0.95;
		% interpolation order : 0 (constant), 1 (linear), 2 (quadratic)
%		order;
		% semivariance model : idw, exponential, Gaussian
		model;
		% maximum distance for the idw method
		% TODO, conflicts with thresh parameter, has similar use
		idw_max_dist = 1000;

		% number of repetitive semivariance estimations
		% TODO is the lsqnonlin residual not enought to access certainty?
		n_var_est_repetitions = 16;
		% number of points used to estiamte the variogram
		n_var_est_samples = 400;
		% initial distance around points to estimate variogram
		% TODO, this distance should be adapted during the estimation process
		var_est_dist = 1000;
		% parameter estimation data
		est;	
		% data of verification run
		verification = struct('n',1000);
	end % properties

	methods

	% constructor
	function obj = Kriging(Rmax,Rmin,order,nverify,aspect_ratio,model)
		% call superclass constructor
		obj = obj@Interpolator(Rmax,Rmin,order,nverify,aspect_ratio);
		obj.model = Model(model);
		obj.order = Order(order);
	end % constructor

	% TODO the existence of this function is due to faulty class design
	% thresh and param and also the cooresponding logic have to go into the model class
	function [dx obj] = ellipsis(obj)
		dx = obj.model.ellipsis(obj.thresh,obj.param);
	end

	end % methods (class Kriging)
	
end % class Kriging

