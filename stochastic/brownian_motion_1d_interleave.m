% Sun 12 Jun 17:22:44 CEST 2022
%
% generate random walk (brownian noise/motion) by interleaving
% TODO clip when number does not match n = 2^p+1
% TODO scale for distance

function a=brownian_noise_interleave(n1,n2,tflag)
	if (nargin()<3)
		tflag = false;
	end

	n = 2.^(nextpow2(n1-1))+1;
	m = n1;	

	% values at end point
	a      = zeros(n1,n2);
	if (~tflag)
		a(n,:) = sqrt(m-1)*randn(1,n2);
	else
		% compute variance
		a(n,:) = (m-1);
	end
	k = 1;
	while (m>2)
		m_ = (m-1)/2+1;
		l = 1;
		r = l+m-1;
		for idx=1:k
			mu = 0.5*(a(l,:) + a(r,:));
			if (~tflag)
				a(l+m_-1,:) =  mu + sqrt((m_-1)/2)*randn(1,n2);
			else
				a(l+m_-1,:) =  mu;
			end
			l = l+m-1;
			r = r+m-1;
		end
		m = m_;
		k = 2*k;
	end
end

