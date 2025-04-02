% Tue 23 Aug 10:49:44 CEST 2022
% Karl Kastner, Berlin
%
%% function [Sc] = lowpass2dpdf_max(f0,p)
%%
%% output: 
%% Sc : maximum of the isotropic two-dimensional lowpass filter
%%      in continuous space
%%
%% input :
%% f0 : frequency scale
%% p  : order
function [Sc] = lowpass2dpdf_max(f0,p)
	Sc = ((2*p-1))./(pi*f0.^2);
end

