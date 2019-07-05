% Mon 11 Jul 10:35:29 CEST 2016
% Karl Kastner, Berlin
%
%% wrap the phase to +/- pi
%
function phi = phasewrap(phi)
	phi = mod(phi+pi,2*pi)-pi;
end

