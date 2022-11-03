% Thu  6 Jan 15:53:21 CET 2022
% Karl Kästner, Berlin
%
%% maximum of the bandpass spectral density
%
function [Sc,IS] = spectral_density_bandpass_continuous_max(fc,p,varargin)
	Sc = spectral_density_bandpass_continuous_scale(fc,p,[],varargin{:});
end

