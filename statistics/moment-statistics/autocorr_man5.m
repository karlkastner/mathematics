% 2015-06-05 14:05:32.609047998 +0200
% Karl Kastner, Berlin
%
%% autocorrellation of the columns of X
%
% function a = autocorr1(X,varargin)
function a = autocorr_man5(X,varargin)
	% a = zeros(n+1,size(X,2));
	a = [];
	for idx=1:size(X,2)
		a(:,idx) = autocorr(X(:,idx),varargin{:});
	end
end

