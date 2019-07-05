% 2016-04-16 21:05:50.029614572 +0200
% Karl Kastner, Berlin
%
%% root mean square value when sample contains nan-values
function x = nanrms(x,varargin)
	x = sqrt(nanmean(x.^2,varargin{:}));
end
