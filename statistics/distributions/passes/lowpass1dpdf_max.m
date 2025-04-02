% Tue 23 Aug 10:49:44 CEST 2022
% Karl Kastner, Berlin
%
%% ouput:
%% Sc : maximum of the density of the spatial, i.e. two-side, low-pass in one dimension
%%
%% input :
%% fc : frequency scale
%% p  : order
function Sc = lowpass1dpdf_mode(f0,p)
	Sc = gamma(2*p)./(f0.*gamma(3/2).*gamma(2*p-1/2));
end

