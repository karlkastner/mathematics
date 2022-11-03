% Fri 10 Sep 16:58:56 CEST 2021
function [mu_S,sd_S,q_S,q_F,pval_r,pval_i] = periodogram_bootstrap(y,m,nb,q,p)
	n  = length(y);
	np = floor(n/m);
	if (nargin()<3 || isempty(nb))
		nb = n-np-1;
		full = false; % true
	else
		full = false;
	end
	w = tukeywin(np);

	S = zeros(n,nb);
	for idx=1:nb
		f = 0;
		% draw m times with replacement
		y__ = zeros(n,1);
		for jdx=1:m
			% sample a random
			if (full)	
				k = idx-1;
			else
				k = randi(n-np+1)-1;
			end
			y_ = zeros(n,1);
			%y_(k+(1:np)) = y(k+(1:np));
			y_j = w.*y(k+(1:np));
			y_j = y_j - mean(y_j);
			y__(k+(1:np)) = y__(k+(1:np)) + y_j;
			%f = f + fft(y_);
		end
%			figure(1)
%			plot(y__); %fft(y_)))
%			pause
		f = fft(y__);
		%rms(f_-f)./rms(f)
		f = f/m;
		F(:,idx) = f;
		S(:,idx) = f.*conj(f);
	end
	% compute the mean
	mu_S = mean(S,2);
	% compute the standard deviation
	sd_S = std(S,[],2);
	% quantiles
	if (nargin()>3)
		q_S   = quantile(S',q)';
		q_F   = quantile(real(F)',q)' + 1i*quantile(imag(F)',q)';
		for idx=1:n
		y    = (1:nb)/(nb+1);
		x    = sort(real(F(idx,:))');
		fdx  = [true; x(2:end)~=x(1:end-1)];
		if (sum(fdx)>1)
		pval_r(idx,1)  = interp1(x(fdx),y(fdx),p,'linear')';
		end
		x    = sort(imag(F(idx,:))');
		fdx  = [true; x(2:end)~=x(1:end-1)];
		if (sum(fdx)>1)
		pval_i(idx,1)  = interp1(x(fdx),y(fdx),p,'linear')';
		end
		end
	end
end
