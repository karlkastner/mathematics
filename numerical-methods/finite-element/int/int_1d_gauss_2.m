% 2012 Apr 27 15:49 (MSK)
% Karl Kästner, Berlin

% coordinates and weights for numerical Gauss quadrature
% 3rd order accurate
function [w, b, flag] = int_1d_gauss_2()

	% quadrature points in baricentric coordinates	
	%a = 0.5*sqrt(1./3) + 0.5
	%b = [   a (1-a);
        %    (1-a)   a];
	% w = [1];
	% b = [0.577350269189626];
	% weights
 	w = [ 0.500000000000000; 0.500000000000000];
	b = [ 0.211324865405187 0.788675134594813;
              0.788675134594813 0.211324865405187];

	% no diagonal mass-matrix
	flag = 0;
end % int_1d_gauss_2()

