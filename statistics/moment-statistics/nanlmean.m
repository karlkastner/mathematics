% Di 7. Apr 11:53:21 CEST 2015
% Karl Kastner, Berlin
%
%% mean of the l-th power of the absolute value of x
function mu = nanlmean(x,l,varargin)
	mu = nanmean(abs(x).^l,varargin{:}).^(1/l);
end
