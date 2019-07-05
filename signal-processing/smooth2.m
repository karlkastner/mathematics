% 2015-06-03 14:52:49.081876531 +0200
% Karl Kastner, Berlin
%% smooth vectos of X
function X = smooth2(X,varargin)
	for idx=1:size(X,2)
		X(:,idx) = smooth(X(:,idx),varargin{:});
	end
end

