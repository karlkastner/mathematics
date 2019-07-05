% Thu 22 Sep 10:46:08 CEST 2016
% Karl Kastner, Berlin
% apply laplacian along columns of x
function x = laplacian(x,n)
% TODO use conv
	for idx=1:n
		x(2:end-1,:) = 0.25*[x(1:end-2,:) + 2*x(2:end-1,:) + x(3:end,:)];
	end
end

