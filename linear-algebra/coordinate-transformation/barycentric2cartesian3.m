% 2016-06-08 11:25:13.371240344 +0200
% Karl Kastner, Berlin
%
%% convert barycentric to cartesian coordinates
% pq : 9xn
function xy = barycentric2cartesian3(a,b,c,pqr)
	xy = zeros(6,size(a,2));
	xy(1,:) = a(1,:).*pqr(1,:) + b(1,:).*pqr(2,:) + c(1,:).*pqr(3,:);
	xy(2,:) = a(2,:).*pqr(1,:) + b(2,:).*pqr(2,:) + c(2,:).*pqr(3,:);
	xy(3,:) = a(1,:).*pqr(4,:) + b(1,:).*pqr(5,:) + c(1,:).*pqr(6,:);
	xy(4,:) = a(2,:).*pqr(4,:) + b(2,:).*pqr(5,:) + c(2,:).*pqr(6,:);
	xy(5,:) = a(1,:).*pqr(7,:) + b(1,:).*pqr(8,:) + c(1,:).*pqr(9,:);
	xy(6,:) = a(2,:).*pqr(7,:) + b(2,:).*pqr(8,:) + c(2,:).*pqr(9,:);
end % pq2xy9

