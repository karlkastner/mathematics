% 2015-08-04 12:27:19.306873735 +0200
% Karl Kastner, Berlin
%
%% impulse response of the ar2 process
%
% function [x A z] = acfar2(s,r1,r2,n)
% function [x A z] = acfar2(r1,r2,n)
%
function [x A z] = acfar2(varargin)
	switch (length(varargin))
	case {3}
		sigma = 1;
		r1 = varargin{1};	
		r2 = varargin{2};	
		n  = varargin{3};
	case {4}
		sigma = varargin{1};
		r1 = varargin{2};	
		r2 = varargin{3};	
		n  = varargin{4};
	otherwise
		error('here')
	end


%	xm2 = 0;
%	xm1 = 1;
	z = roots([1 -r1 -r2])
	%z = 1./z
	% syms z1 z2 r1
	% why r1/(1-r2) and not simply r1 ?
	A = [    1 1;
             z(1), z(2)] \ [1; r1/(1-r2)]
	% A = [  (r1 - z2)/(z1 - z2), -(r1 - z1)/(z1 - z2)]
	%A    = [(r1 - z(2))/(z(1) - z(2)); -(r1 - z(1))/(z(1) - z(2))]
	x    = zeros(n,1);
	N = (0:n-1)';
	x = A(1)*z(1).^N + A(2)*z(2).^N;
	x = sigma*x;
%	x(:,2) = A(1)*z(1).^N;
%	x(1) = 1;
%	x(2,1) = r1/(1-r2);
%	for idx = 3:n
%		x(idx,1) = r1*x(idx-1) + r2*x(idx-2); % + c*0;
%	end
%	close all
%	norm([x-x_])
%	plot([x x_]);
end

