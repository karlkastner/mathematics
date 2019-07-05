% Di 10. Nov 17:16:14 CET 2015
% Karl Kastner, Berlin
%
%% interpolation in streamwise coordinates
%
% This gives similar result to setting aspect ratio for N to infinity,
% but not quite,as the input point set is not dense (scale for N to infinity does not work)
%
% input:
%
%	sS  : S-coordinate of source points
%       tN  : N-coordinate of source points
%	tV  : value at source points
%       sdx : segment index of source points
%       tS  : S-coordinate of target points
%	tN  : N-coordinate of target points
%	tdx : segment index pf target points
%	
% output:
%	tV  : value interpolated to target points

% TODO at bifurcations: follow deeper channel, or follow channel that minimises difference in cross section
% TODO special treatment to outside points
%      add artifical lines from nmax/w to 0.5 and nmin/w to -0.5
function tV = interp_sn2(sS,sN,sdx,sV,tS,tN,tdx,S_range,L_max,order)

	% build an source point index for all segments
	nseg = max(max(sdx),max(tdx));
	sid  = cell(nseg,1);
	for idx=1:nseg
		sid{idx} = find(idx == sdx);
	end

	% allocate memory
	ns     = length(sS);
	nt     = length(tS);
	tV     = NaN(nt,size(sV,2));
	L_max2 = L_max*L_max;
	
	% for each target point
	for idx=1:nt
		% get all points within range
		% TODO, only in the same segment
		if (tdx(idx) <= 0)
			continue;
		end

		sdx = sid{tdx(idx)};
		%[tV(idx,:)] = interp_sn_();
		tV = interp_sn_(idx,sdx,sS,sN,tS,sV,tN,tV,ns,order,L_max2);
	end % for idx

end % function

