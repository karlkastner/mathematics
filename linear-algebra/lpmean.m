% 2015-03-03 14:12:53.340974566 +0100
% Karl Kastner, Berlin
%
%% mean of pth-power of a
function n = lpmean(A,p,dim)
	if (nargin()<2)
		p = 2;
	end
	if (nargin()<3)
		dim = 1;
	end
	n = (mean(abs(A).^p,dim)).^(1/p);
end

