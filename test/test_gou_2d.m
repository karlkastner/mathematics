% Tue 28 Nov 10:54:33 CET 2023

% variance seems not to be correct, even if order is increased to keep sampling the same
% corr is too large so var is too low
% same proble in 1d and 2d
% works when correlation is zero

geo = false;

	L     = 0.125*1024/4;
	n     = 1024/4;
	dynamic_order = true;
	order = 10;
	theta = 0.35;
	pmu   = 0.2;
	psd   = 0.1;
	idx   = 0;	
	mu    = [];
	sd    = [];
	sk    = [];
	r     = [];
	 while (n>16)
	 idx=idx+1;
	 disp(idx)
	if (geo)
		 y = geometric_ar1_2d_grid_cell_averaged_generate(lmu,lsd,theta,L*[1,1],n*[1,1],order);
	else
		 y = ar1_2d_grid_cell_averaged_generate(pmu,psd,theta,L*[1,1],n*[1,1],order);
	end
	if (1 == idx)
		x = y;
	end

	 N(idx) = n;
	 mu(1,idx)  = mean(x(:));
	 sd(1,idx) = std(x(:));
	 sk(1,idx) = skewness(x(:));
	 r(1,idx)  = corr(cvec(x(1:end-1)),cvec(x(2:end)));
	 mu(2,idx) = mean(y(:));
	 sd(2,idx) = std(y(:));
	 sk(2,idx) = skewness(y(:));
	 r(2,idx)  = corr(cvec(y(1:end-1)),cvec(y(2:end)));
if (0)
	figure(1);
	subplot(2,5,idx)
	R = ifft2(abs(fft2(x-mean(x(:)))).^2);
	R = R/R(1);
	R_ = ifft2(abs(fft2(y-mean(y(:)))).^2);
	R_ = R_/R_(1);
	plot([R(:,1),R_(:,1)])
end

	 A = downsampling_matrix(N(idx),'pairwise');
	 x = (A'*x*A);
	 n=n/2;
	if (dynamic_order)
	 order = 2*order;
	end
	 end
mean(x)
L./N
mu
sd
sk
r
sd(:,2:end)./sd(:,1:end-1) 

