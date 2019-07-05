% Mi 8. Apr 16:07:28 CEST 2015
% Karl Kastner, Berlin
%
%% amplitude and standard deviation of the amplitude of a frequency component
%% represented by a peak in the fourier domain
%
%% input :
%% h : peak height
%% w : peak width at half height
%%
%% output:
%% a : amplitude in real space
%% s : standard deviation of the frequency (!)
%
% TODO is this standard deviation of standard error?
function [a, s] = apmlitude_from_peak(h,w)
	% standard deviation from peak width
	s = w.^2/(8*log(2));
	% amplitude from peak width
	a = h./normpdf(0,0,s);	
end

