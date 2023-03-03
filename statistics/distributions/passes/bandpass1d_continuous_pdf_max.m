% Thu  6 Jan 15:53:21 CET 2022
% Karl Kästner, Berlin
%
%% maximum of the bandpass spectral density
%
function [Sc,IS] = bandpass1d_continuous_pdf_max(fc,p,varargin)
	Sc = bandpass1d_continuous_pdf_scale(fc,p,[],varargin{:});
end

