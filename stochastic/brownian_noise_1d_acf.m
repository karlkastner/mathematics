% Wed 26 Jan 17:14:34 CET 2022
%
% generate brownian noise from its autocorrelation
%
% note: this is from the definition of fractional brownian noise
%       and not computationally efficient
% c.f. wikipedia: Fractional Brownian motion
% Dietrich and newham: how to use, brownian is not stationary?
function e = brownian_noise_1d_acf(L,n,e0)
	if (nargin()<3)
		if (length(n)<2)
			n(2) = 1;
		end
		e0 = randn(n);
	end	
	t = L*(0:n(1)-1)'/n(1);
	t = innerspace(0,L,n(1))';
	%t = L*(1:n)'/n;
	% define the covariance matrix
	% = min(t,s)
	s = t';
	R = 0.5*(s + t - abs(t-s));
	S = sqrtm(R);
%[v,e] = eig(R);
%R
%[flipud(sort(diag(e))) , (1./(1:n(1))'.^2)]
%d = diag(R);
%d=[d;flipud(d(2:end))]
%fft(d)
%pause
	% gaussian noise
	e  = S*e0;
end


