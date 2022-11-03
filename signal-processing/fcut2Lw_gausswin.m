% 2022-05-17 11:01:33.529366562 +0200
% cutoff frequency to window length
function Lw = fcut2Lw_gausswin(fcut)
	% Lw = sqrt(-log(0.5)/4)*1./(pi*fcut);
	% this is it's own inverse
	Lw = fcut_gausswin(fcut);
end

