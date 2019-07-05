% 2014-12-02 18:43:34.593549779 +0100

%% limit to +- pi
% TODO use wrapphase
function phase = pilim(phase)
	%phase = phase + 2*pi*(phase < -pi) - 2*pi*(phase > pi);
	phase = phase + 2*pi*(phase < 0) - 2*pi*(phase > 2*pi);
end

