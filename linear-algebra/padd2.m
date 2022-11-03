% Fri 10 Dec 11:25:14 CET 2021
% Karl Kästner, Berlin
%
%% padd values around a 2d (image) matrix, constant exprapolation
%
% function x = padd2(x,nf)
function x = padd2(x,nf,z)
	if (nargin()<3)
		z = 1;
	end
	if (length(nf)==1)
		nf(2) = nf(1);
	end
	if (z)
	x = [repmat(x(1,:),nf(1),1); x; repmat(x(end,:),nf(1),1)];
	x = [repmat(x(:,1),1,nf(2)),x,repmat(x(:,end),1,nf(2))];
	else
	x = [0.*repmat(x(1,:),nf(1),1); x; 0.*repmat(x(end,:),nf(1),1)];
	x = [0.*repmat(x(:,1),1,nf(2)),x, 0.*repmat(x(:,end),1,nf(2))];
	end
end

