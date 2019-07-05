% Mo 6. Apr 15:10:09 CEST 2015
% Karl Kastner, Berlin
%
%% mean of x.^l, not of abs
function mu = lmean(x,l,varargin)
	mu = mean(x.^l,varargin{:}).^(1/l);
end
