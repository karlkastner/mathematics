% Mon 11 Jul 09:45:24 CEST 2016
% 2015-05-16 14:54:53.828308851 +0200
% Karl Kastner, Berlin
%
%% fill nan-values in vector with gaps
% function y = fixnan(x,y,dx_max,varargin)
%
function y = fixnan(x,y,dx_max,varargin)
	y = interp1_limited(x,y,x,dx_max,varargin{:});
%	x = cvec(x);
%	if (isvector(y))
%		y = cvec(y);
%	end
%	if (length(x) ~= size(y,1))
%		error('Length of x has to match number of rows in y');
%	end
%
%	% centre x to reduce cancellation error
%	% (especially relevant for single precission and cubic interpolation)
%	xc = 0.5*(x(1)+x(end));
%	x  = x-xc;
%
%	for idx=1:size(y,2)
%		% interpolate values over gaps	
%		flag = isnan(y(:,idx));
%		if (sum(~flag)>1)
%			% interpolate
%			y(flag,idx) = interp1(x(~flag),y(~flag,idx),x(flag),varargin{:});
%
%			% get maximum difference along x
%			dx       = zeros(size(y,1),1);
%			dx(flag) = interp1(x(~flag),x(~flag),x(flag),'nearest')-x(flag);
%	
%			% invalidate sample in gapes more than dx_max apart
%			%flag = (abs(t_nearest-x) > dx_max);
%			y(abs(dx)>dx_max,idx) = NaN;
%		end % if
%	end % for idx
end % fixnan

