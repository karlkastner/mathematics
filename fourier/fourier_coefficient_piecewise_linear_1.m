% Tue 15 May 11:56:36 CEST 2018
% Karl Kastner, Berlin
%
%% fourier series coefficients of a piecewise linear function
%% (not coefficient of discrete fourier transform)
%% function can be discontinuous between intervals
%% scales domain length to 2pi
%%
%% input :
%% X : end points of piecewise linear function
%% Y : values at end points
%% 
%% output :
%% ab : coefficients for frequency components
%
function ab = fourier_coefficient_piecewise_linear(X,Y,kk)
	L = X(end);
	X = 2*pi*X/L;

%ab = zeros(2,length(kk));
	for jd=1:length(kk);
	k   = kk(jd);
	A   = 0;
	rhs = 0;
	for id=1:length(X)-1
		% right - left
		% TODO, this cancel for all but the first and last
		x0 = X(id);
		x1 = X(id+1);
		y0 = Y(1,id);
		y1 = Y(2,id);
		dy = (y1-y0)./(x1-x0);
		y0 = y0 - dy*x0;
		if (0 == k)
		A   = A + ( [1 0; 0 2*(x1-x0)] );
		rhs = rhs + [0; 
			     (x1*(6*y0 + 3*dy*x1))/3  - (x0*(6*y0 + 3*dy*x0))/3 ];
		else
		A = A - ( [ -(3*sin(2*k*x0) - 6*k*x0)*k/6,                 -cos(k*x0)^2*k;
		                             -cos(k*x0)^2*k, (3*sin(2*k*x0) + 6*k*x0)*k/6] ...
		        - [ -(3*sin(2*k*x1) - 6*k*x1)*k/6,                 -k*cos(k*x1)^2;
		                             -cos(k*x1)^2*k, (3*sin(2*k*x1) + 6*k*x1)*k/6] );
	
		rhs = rhs - ( [  2*dy*sin(k*x0) + -2*(y0 + x0*dy)*k*cos(k*x0);
			         2*k*(y0 + x0*dy)*sin(k*x0) + 2*dy*cos(k*x0) ] ...
			    - [  2*dy*sin(k*x1) + -2*(y0 + x1*dy)*k*cos(k*x1);
			         2*k*(y0 + x1*dy)*sin(k*x1) + 2*dy*cos(k*x1) ] );
		end
	end % for id
	ab(:,jd) = A \ rhs;
	end % for jd
end % pwl2fourier

