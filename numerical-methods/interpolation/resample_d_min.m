% Mon 12 Dec 20:23:00 CET 2016
%% resample a function
function [xr varargout] = resample_d_min(x,dx_min,varargin)
	last = zeros(size(x));
	last(1) = 1;
	n = 1;
	tic()
	x = double(x);

	fdx = isfinite(x);	
	if (any(~fdx))
		warning('non-infinite sample in input vector ignored');
		x = x(fdx);
	end

	dx_min = double(dx_min(1));
	id = javaMethod('resample_d_min','Resample',x,dx_min);
	toc()

	xr = accumarray(id,x,[id(end) 1],@mean);

	% TODO, do this in java
%	for idx=2:length(x)-1
%		if (abs(x(idx)-x(last(n)))>=dx_min)
%			%idx/length(x)
%			n = n+1;
%			last(n) = idx;
%		end
%	end
%	% tail rule
%	if (abs(x(end)-x(last(n)))<0.5*d_min)
%		last(n) = length(x);
%	else
%		n=n+1;
%		last(n) = length(x);
%	end

	% resample x
%	xr = zeros(n-1,1);
%	for idx=1:n-1
%		xr(idx) = mean(x(last(idx):last(idx+1)));
%	end
	varargout = {};
	% resample other
	for jdx=1:length(varargin)
		v  = varargin{jdx};
		v  = v(fdx);
	%	vr = zeros(n-1,1);
		vr = accumarray(id,v,[id(end) 1],@mean);
%		for idx=1:n-1
%			vr(idx) = mean(v(last(idx):last(idx+1)));
%		end
		varargout{jdx} = vr;
	end
end

