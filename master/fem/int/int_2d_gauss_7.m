% Sun May 20 20:45:24 MSK 2012
% Karl Kästner, Berlin

% weigts for 4th order accurate gauss quadrature
function [w p flag] = int_2d_gauss_7()
%	warning('does not work in combinatiion with 15-point scheme');	

	% integration weights
	w = [	0.225000000000000;
		0.125939180544827; 0.125939180544827; 0.125939180544827;
		0.132394152788506; 0.132394152788506; 0.132394152788506 ];

	% integration points in barycentric coordiantes
	p = [	0.333333333333333 0.333333333333333 0.333333333333333
		0.797426985353087 0.101286507323456 0.101286507323456
		0.101286507323456 0.797426985353087 0.101286507323456
		0.101286507323456 0.101286507323456 0.797426985353087
		0.059715871789770 0.470142064105115 0.470142064105115
		0.470142064105115 0.059715871789770 0.470142064105115
		0.470142064105115 0.470142064105115 0.059715871789770 ];

	% mass matrix will not be diagonal
	flag = 0;
end % int_2d_gauss_7

