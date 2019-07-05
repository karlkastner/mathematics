% Wed 15 Jun 15:26:51 CEST 2016
%% smooth subsequent values along a curve such that
%%	v(x0+dx) < v(x0) + (ratio-1)*dx
%% if v is the edge length in a resampled polygon, then v_i/v_(i+1) < ratio
%% 	ratio^1 = exp(a*1)
function val = limit_by_distance_1d(X,Y,val,ratio)
	X   = cvec(X);
	Y   = cvec(Y);
	val = cvec(val);
	
	% split in segments
	fdx = [1;cvec(find(isnan(X)));length(X)];
	seg = [fdx(1:end-1)+1,fdx(2:end)-1];
	% determine segment length
	dS  = NaN(size(X));
	for idx=1:size(seg,1)
		fdx = seg(idx,1):seg(idx,2);
		if (~isempty(fdx))
			% TODO use left and right distance, not cdiff
			dS_     = hypot(cdiff([X(fdx(end));X(fdx);X(fdx(1))]),cdiff([Y(fdx(end));Y(fdx);Y(fdx(1))]));
			dS_     = dS_(2:end-1);
			dS(fdx) = dS_;
		end % if
	end % for idx

	while (true)
		val_old  = val;	
		for idx=1:1:size(seg,1)
			fdx = seg(idx,1):seg(idx,2);
			if (~isempty(fdx))
				% right
				val(fdx)  = min(val(fdx),([val(fdx(2:end));val(fdx(1))]+(ratio-1)*dS(fdx)));
				% left
				val(fdx)  = min(val(fdx),([val(fdx(end));val(fdx(1:end-1))]+(ratio-1)*dS(fdx)));
			end % if
		end % for idx
		if ( 0 == max(abs(val-val_old)))
			break;
		end
	end % while true
end % limit_by_distance_1d

