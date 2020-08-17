
	% resampling of solution
	if (isfield(opt,'xr') && ~isempty(opt.xr))
		xr = opt.xr;
		k  = 1;
		yr = zeros(size(xr));
		for idx=1:length(xr)
			while(xr(idx) > x(k+1) && k+1<nx)
				k = k+1;
			end % while
			% expand
			dx   = xr(idx) - xc(k);
			yr(idx) =   ypm(3*(k-1)+1)*myexp(dx*l(k,1)) ...
		                  + ypm(3*(k-1)+3)*myexp(dx*l(k,2));
		end % for idx
	end % if opt.xr

