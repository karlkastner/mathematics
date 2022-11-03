% 2015-08-11 16:49:43.834138528 +0200
% Karl Kastner, Berlin
%% moving median filter with special treatment of boundaries

function Y = medfilt1_man(X,nf)
	if (isvector(X))
		X = cvec(X);
	end
	% delta
	d = fix((nf-1)/2);
	
	n = size(X,1);
	Y = zeros(size(X));
	for idx = 1:d
		% ramp up
		Y(idx,:) = median(X(1:2*idx-1,:));
		% ramp down
		Y(n-idx+1,:) = median(X(n-2*idx+2:end,:));
	end
	for idx=d+1:n-d
		Y(idx,:) = median(X(idx-d:idx+d,:));
	end
end
