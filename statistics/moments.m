% 2024-05-14 21:49:11.816930914 +0200
% Karl Kastner, Berlin
% non-central moment
function mu   = moments(x,n,iscentral)
	mu    = zeros(n,1);
	mu(1) = mean(x);
	if (iscentral)
	x     = x - mu(1);
	end
	xk    = x; 
	for k=2:n
		xk = x.*xk;
		mu(k) = mean(xk);
	end	
end

