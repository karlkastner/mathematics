function x = bandpass2d_iso(x,p,k)
	p = p;
%	p = 0.4;
	for idx=1:k
	x = lp2d_iso(x,p,k);
	%x = x-lp2d_iso(x,p,k);
	end
end

function x = lp2d_iso(x,p,k)
	L = sqrt(2);
	L = 1;
	wy = 1;
	wx = 1;
	wxy = 0;
p=0.8;
	for jdx=1:2
	for idx=2:size(x,1)
		x(idx,:) = (1-p)*x(idx,:) + p*x(idx-1,:);
	end
	for idx=2:size(x,2)
		x(:,idx) = (1-p)*x(:,idx) + p*x(:,idx-1);
	end
	for idx=2:size(x,1)
		x(idx,2:end) = (1-(p.^sqrt(2)-p^2))*x(idx,2:end) + (p.^sqrt(2)-p^2)*x(idx-1,1:end-1);
		%x(idx,2:end) = (1-p)*x(idx,2:end) + (p)*x(idx-1,1:end-1);
	end
		x = fliplr(flipud(x));
	end

if (0)
	for idx=2:size(x,1)
         for jdx=2:size(x,2)
		x(idx,jdx) = (1-(wx+wy)*p-(wxy*p.^L))*x(idx,jdx)  ...
				+ wx*p*x(idx-1,jdx) ...
				+ wy*p*x(idx,jdx-1) ...
				+ (wxy*p.^L)*x(idx-1,jdx-1);
	 end
	end
	for idx=2:size(x,1)
         for jdx=size(x,2)-1:-1:1
		x(idx,jdx) = (1-(wx+wy)*p-(wxy*p.^L))*x(idx,jdx)  ...
				+ wx*p*x(idx-1,jdx) ...
				+ wy*p*x(idx,jdx+1) ...
				+ (wxy*p.^L)*x(idx-1,jdx+1);
	 end
	end
end
end
