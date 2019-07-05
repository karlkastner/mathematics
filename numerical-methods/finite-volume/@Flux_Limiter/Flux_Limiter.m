% Sun 12 Nov 18:33:42 CET 2017
%% class of flux limiters
classdef Flux_Limiter % not inherited from handle
	methods (Static)
		phi = beam_warming(theta);
		phi = fromm(theta);
		phi = lax_wendroff(theta);
		phi = minmod(theta);
		phi = monotized_central(theta);
		phi = superbee2(theta);
		phi = superbee(theta);
		phi = upwind(theta);
		phi = vanLeer(theta);
	end % static methods
	methods
		function obj = Flux_Limiter()
		end
	end
end % class Flux_Limiter

