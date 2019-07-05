% Di 26. Jan 14:03:32 CET 2016
% Karl Kastner, Berlin
%
%% batch interpolation by the fourier interpolation
function [tF tErr] = interp_fourier_batch(sS, sNrel, sdx, sF, tS, tNrel, tdx)

	% allocate memory
	%tF = NaN(numel(tS),size(sF,2),class(sF));
	tF = NaN(size(tS),class(sF));
	tErr = NaN(size(tS),class(sF));

	% list of segments to be processed
	seg = unique(tdx);

%	sid = cell(length(seg),1);
	% interpolate points for each segment individually
	for idx=rvec(seg)
		sid = find(idx == sdx);
		tid = find(idx == tdx(:,1));
		if (~isempty(sid) && ~isempty(tid))
			[tF(tid,:) serri] = interp_fourier(sS(sid),sNrel(sid),sF(sid),tS(tid,:),tNrel(tid,:));
			tErr(tid,:) = serri;
		end
	end
end

