% Fri 10 Dec 11:25:14 CET 2021
% Karl Kästner, Berlin
%
%% padd values around a 2d (image) matrix, constant exprapolation
%
% function x = padd2(x,nf)
function x = padd2(x,nf)
	x = [repmat(x(:,1),1,nf),x,repmat(x(:,end),1,nf)];
	x = [repmat(x(1,:),nf,1); x; repmat(x(end,:),nf,1)];
end
