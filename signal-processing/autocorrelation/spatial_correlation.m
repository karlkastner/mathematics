% 2024-05-18 17:46:12.353864484 +0200
% Karl Kastner, Berlin
%% spatial autocorrelation, determined from the first lag
%% this is only appropriate, when the acf is first order autocorrelated
%% without moving average of higher order ar-terms
%% it would be better to determine it by the first distance where the correlation
%% drops to exp(-1)
% TODO compute distance L
function [cx,cy] = spatial_correlation(S)
	s  = size(S);
	S_ = [S(end,:);S(1:end-1,:)];
	cx = nancorr(S(:),S_(:));
	cx = cx(1,2);
	if (nargout()>1)
		S = S';
		S_ = [S(end,:);S(1:end-1,:)];
		cy = nancorr(S(:),S_(:));
		cy = cy(1,2);
	end
end

