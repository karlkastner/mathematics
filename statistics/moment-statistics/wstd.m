% So 19. Jul 12:45:11 CEST 2015
% Karl Kastner, Berlin
%
%% weighed standard deviation
%
% function [sd] = wstd(w,x)
function [sd] = wstd(w,x,varargin)
	sd = sqrt(wvar(w,x,varargin{:}));
end % wstd

