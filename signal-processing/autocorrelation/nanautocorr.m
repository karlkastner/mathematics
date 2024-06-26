% 2015-04-09 19:41:16.007183956 +0200
% Karl Kastner, Berlin
%
%% autocorrelation with nan-values
% TODO : make frequencies orthogonal
%
function out = nanautocorr(data,nlag)
	if (nargin() < 2)
		nlag = 20;	
	end
	nlag = min(nlag,length(data)-1);
	out=zeros(nlag,1);
	out(1)=1;
	for idx=2:nlag
	    out(idx)=corr(data(idx:end),data(1:end-idx+1),'rows','complete');
	end
end

