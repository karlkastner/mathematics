% Sun 12 Jun 17:22:44 CEST 2022
%
% generate brownian surface by interleaving
% aka Gaussian free field, Brownian sheet
%
% TODO clip when number does not match n = 2^p+1
% TODO scale for L ~= n
% TODO non-square domains

function a=brownian_noise_interleave(n1,test)
	if (nargin()<2)
		test = false;
	end

	n = 2.^(nextpow2(n1-1))+1;
	m = n1;	
	n2 = n1;

	% values at corners
	a      = zeros(n1,n1);
	if (~test)
		a(n,1) = sqrt(m-1)*randn();
		a(1,n) = sqrt(m-1)*randn();
		% TODO, does the corner depend on the others?
		%a(n,n) = sqrt(2*(m-1))*randn();
		a(n,n) = sqrt((m-1))*randn();
	else
		a(n,1) = 1;
		a(1,n) = 1;
		a(n,n) = 1; %sqrt(2);
	end

	k = 1;
	while (m>2)
		m_ = (m-1)/2+1;
		% diagonal interleaving
		l = 1;
		r = l+m-1;
		for idx=1:k
		 b = 1;
		 t = b+m-1;
		 for jdx=1:k
			mu = 0.25*(a(l,t) + a(r,t) + a(l,b) + a(r,b));
			if (~test)
				a(l+m_-1,b+m_-1) =  mu + sqrt(2*(m_-1)/2)*randn();
			else
				a(l+m_-1,b+m_-1) =  mu;
			end
			t = t+m-1;
			b = b+m-1;
		 end
		 l = l+m-1;
		 r = r+m-1;
		end
if (0)
		% vertical interleaving
		i = m_;
		for idx=1:k
		 j = m_;
		 for jdx=1:k
			%	    l,r,t,b
			mu = 0.25*(a(i,j) + a(i,j) + a(i,j) + a(i,j));
			a(l+m_-1,b) =  mu + sqrt((m_-1)/2)*randn();
			t = t+m-1;
			b = b+m-1;
		 end
		 l = l+m-1;
		 r = r+m-1;
		end
		% horizontal interleaving
a
pause
end
		m = m_;
		k = 2*k;
	end
end

