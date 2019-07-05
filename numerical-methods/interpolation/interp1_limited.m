% Sa 20. Feb 10:44:47 CET 2016
% Karl Kastner, Berlin
%% interpolate values, but not beyond a certain distance
%% this function is idempotent, i.e. it will not extrapolate over into gaps
%% exceedint the limit and thus not spuriously extend the series when called a second time on the same data
function yi = interp1_limited(x,y,xi,dx_max,extrapflag,varargin)
	if (isempty(dx_max))
		dx_max = inf;
	end

	x  = cvec(x);
	xi = cvec(xi);
	if (isvector(y))
		y = cvec(y);
	end
	if (length(x) ~= size(y,1))
		error('Length of x has to match number of rows in y');
	end
	% make unique and sort
	[x, sdx] = unique(x);

	% centre x to reduce cancellation error
	% (especially relevant for single precission and cubic interpolation)
	xc = 0.5*(x(1)+x(end));
	x  = x  - xc; 
	xi = xi - xc; 

	y       = y(sdx,:);
	% allocate memory
	yi = NaN(length(xi),size(y,2));
	if (length(x) > 1)
		for idx=1:size(y,2)
			fdx       = isfinite(y(:,idx));
			% find is required, as later indices of fdx are accessesed
			fdx       = find(fdx);
			nx        = length(fdx);
			if (nx>1)
				id   = 1:nx;

				% interpolate the value
				yi(:,idx) = interp1(x(fdx),y(fdx,idx),xi,varargin{:});

				% interpolate the index
				idi  = interp1(x(fdx),id,xi,'linear','extrap');
				% note: idil == idir in case of exact match
				idil = floor(idi);
				idir = ceil(idi);
				% no extrapolation
				% extrapolation is dangerous if this function is applied recursively
				fdxi =   idil > 0 & idil <= nx ...
                                       & idir > 0 & idir <= nx;
				%dx_ = NaN(size(x)); dx_(fdxi) = xi(fdxi) - x(fdx(idil(fdxi)));
				%fdxi(fdxi) =   abs(xi(fdxi) - x(fdx(idil(fdxi)))) < dx_max ...
				%             & abs(xi(fdxi) - x(fdx(idir(fdxi)))) < dx_max;
				fdxi(fdxi) =   min(abs(xi(fdxi) - x(fdx(idil(fdxi)))), ...
                                                   abs(xi(fdxi) - x(fdx(idir(fdxi))))) <= dx_max/2;
				if (nargin() > 5 && ~isnan(extrapflag))
					% TODO quick fix for extrapolation
					fdxi(1) = true;
					fdxi(end) = true;
				end
				yi(~fdxi) = NaN;
			end
		end % for
	end % if
end % interp1_limited

