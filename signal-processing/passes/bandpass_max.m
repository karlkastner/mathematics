% 2021-07-02 20:24:42.600181396 +0200
% location, where bandpass reaches maximum
function bp_max = bandpass_max(f0,fmax)
	bp_max = fmax*(1 - 1/pi*acos((f0.^2 - 4*f0*fmax + fmax.^2)/(2*f0*fmax)));
end

