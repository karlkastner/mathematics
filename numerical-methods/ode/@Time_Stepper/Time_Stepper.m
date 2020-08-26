% Thu 20 Aug 10:50:22 +08 2020
classdef Time_Stepper < handle
	properties
		% start and end time of similartion (seconds) 
		Ti


		% numerical scheme for time stepping
		scheme = 'upwind';
		
		% number of terms for multi-step and runge-kutta
		order = 1;

		% factor for courant-friedrich-lewis stability condition
		cfl = 0.99;

		% interval at which the bed level is written to the output (seconds)
		% if 0, write always
		dto = 0;

		% stopping criteria for implicit schemes
		reltol = sqrt(eps);
	
		n_realloc = 100;
	end
	methods

	end % methods
end
