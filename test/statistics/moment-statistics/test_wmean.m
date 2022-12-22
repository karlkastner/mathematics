% 2018-06-05 09:40:11.906487406 +0200
function test_wmean()

	w = ones(1,10);
	x = 1:10;
	
	mu1 = wmean(w,x);

	res(1) = mu - 5.5;

	mu2 = wmean(x,w);
	
	res(2) = mu - 1

	n = 10;
	x = 1:n;
	w = randi(10,1,n);
	x_ = [];
	for idx=1:n
		x_  = [x_,repmat(idx,1,k(idx))];
	end
	mu3 = wmean(w,x)
	mu3_ = mean(x_)
	res(3) = mu3-mu3_;

	fail = max(abs(res))>sqrt(eps)

%	n=1e2;
% 	x = [ones(n-1,1); n];
%	w=ones(n,1);
%	for idx=1:10
%		mu=wmean(w,x)
%		s=(x-mu).^2;
%		w = 1./(mean(s)+s);
%	end
end

