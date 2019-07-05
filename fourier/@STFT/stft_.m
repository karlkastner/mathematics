% Fri  5 May 13:32:53 CEST 2017
% Karl Kastner, Berlin
%% static wrapper for STFT
function [tc, c, amp, phase, serr, s] = stft(val,dt,Ti,T,winfun)
	s = STFT('T',T,'Ti',Ti,'winfun',winfun);
	s.transform(val,dt);
	tc    = s.tc;
	c     = s.coeff;
	amp   = s.amplitude;
	phase = s.phase;
	serr  = NaN; %s.serr;
end

