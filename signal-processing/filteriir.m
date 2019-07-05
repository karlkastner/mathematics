% Thu 31 May 09:45:03 CEST 2018
%
%% filter adcp t-n data over time
%%
%% v : nz,nt   : values to be filtered
%% H : nt,1    : depth of ensemble
%% last : nt,1 : last bin above bottom that can be sampled without side lobe interference
%% nf : scalar : number of reweighted iterations
%%
%% when samples 
%% - distance to bed is reference (advantageous for near-bed suspended transport)
%% TODO for wash load: distance to surface is more relevant
%% interpolate depending on z
%%
%% when depth changes, neighbouring indices do not correspond to same relative position in the water column
%% relative poisition in the colum (s-coordinate) smoothes values
%% near the bed: absolute distance to bed is chosen
%% near surface: absolute distance to surface is chosen
%% -> cubic transformation of index
%%
%% faster and avoid alising (smoothing along z)
%%	resample ensemble to same number of bins in S -> filter -> resample back
%%	use nonlinear transform z-s coordinates
%% -> resampling has to be local (Hi -> H-filtered)
%%
%% filtered profile coordinates to sample coordinates
%% 	zf -> zi (special transform)
%% corresponding indices and fractions
%% filtration step (update of hf and vf)
%% sample coordinates to updated profile coordinates
%% (the inverse step is actually not necessary)
%% write filtered value
%
function [vf,dvf,q] = filteriir(v,H,last,p,nf)
	nens = size(v,2);
	nbin = size(v,1);
	id   = (1:nbin)';

	% apply rangemask
	msk     = isfinite(v) & bsxfun(@le,id,last');
	v(~msk) = 0;

	s_bar  = 1;
	s      = ones(size(v),class(v));
	vl     = NaN(size(v),class(v));
	vr     = NaN(size(v),class(v));

	% iterate for reweighted filtering
	for jdx=1:nf
		vl = filter(v,vl,+1);
		vr = filter(v,vr,-1);
	
		% average
		% TODO filter from right to left
		vf = 0.5*(vl + vr);
	
		% estimate error variance of the unfiltered data
		if (nf > 1)
			s     = abs(v-vf);
			s_bar = sqrt(nanmedian(s(:).^2));
		end
	end % for jdx

	% error variance of the filtered data
	if (nargout() > 1)
		%dvf = vl-vr;
		%dvf = vl-vr;
		dvf = vl(:,1:end-1) - vr(:,2:end);
	end

function [id1,id2,ql] = build_index(H,last)
	H    = rvec(H);
	last = rvec(last);
	% matrix for interpolation from left to right
	% corresponding bin index in previous ensemble
	idl  = id*([H(1);H(1:end-1)]./H);
	% integer part
	id1  = floor(idl);
	% fractional part
	ql    = idl-id1;
	% limit
	id1 = bsxfun(@min,max(1,id1),([last(1);last(1:end-1)]));
	id2 = bsxfun(@min,id1+1,([last(1);last(1:end-1)]));
end

% filter from left to right
function vf = filter(v,vf,d)
	if (d < 0)
		r = [1,1,nens];
	else
		r = [nens,-1,1];
	end
	for idx=r(1):r(2):r(3)
		vlast_ = [];
		if (1==idx || nens == idx || 0 == last(idx+d))
			% initialise
%			vlast = v(:,idx);
			vlast = zeros(nbin,1);
			%fdx   = isnan(vlast);
			%vlast(fdx) = mean(vlast(~fdx));
		else

		% index of last
		mode = 'zb';		
		switch (mode)
		case {'S'}
			idl  = id*(H(idx+d)./H(idx));
			% integer part
			id1  = floor(idl);
			% fractional part
			q    = idl-id1;
			% limit
			id1 = min(max(1,id1),last(idx+d));
			id2 = min(id1+1,last(idx+d));
			% interpolate
			%vlast = (1-q(:,idx)).*vlast(id1(:,idx)) ...
			%	 + q(:,idx).*vlast(id2(:,idx));
			vlast = (1-q).*vlast(id1) ...
				  + q.*vlast(id2);
		case {'zb'}
		        % keep-near bed profile the same
			% note that if there is a shift along the vertical,
			% an attenuation factor has to be (un)applied
			%delta = 0.97;
			delta = nanmean(diff(vlast(1:last(idx+d))));
			if (isnan(delta)) delta = 0; end
			dd = (last(idx+d)-last(idx));
			% id1 = id + last(idx+d)-last(idx);			
			id1 = id + dd;			
			id1 = max(1,id1);
			id1 = min(id1,last(idx+d));
			vlast = vlast(id1);
			vlast = vlast-delta*dd;
		case {'zs'}
		case {'cubic'}
		otherwise
		end % % switch mode

		% specify weight depending on the variance (reweighted)
		% and low pass weight
		% invalid samples have zero weight
		w = p .* msk(:,idx).* 2.*s_bar./(s_bar + s(:,idx));

		% filter
		vlast = w.*v(:,idx) + (1-w).*vlast;
		if (last(idx)>0)
			vlast(last(idx)+1:end) = NaN;
		end
		end

		% store
		vf(1:last(idx),idx) = vlast(1:last(idx));
	end % for idx
end % filter

end % filter iir

