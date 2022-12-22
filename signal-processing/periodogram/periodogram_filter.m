% 2022-05-17 16:07:12.676212489 +0200
% 2022-05-17 16:06:17.459013111 +0200
% estimate the spectral density by filtering (smoothing) the periodogram
% spectral density estimated with a guassian window
function S = periodogram_filter(b,L,m,type)
	if (isvector(b))
		b = cvec(b);
	end
	n = size(b,1);
	switch (type)
	case {'gauss'}
		Lw_spectral = m/L;
		Lw_gauss = Lw_spectral*fcut_gausswin(1)/rectwin_cutoff_frequency(1);
		Ls  = n/L;
		Sf  = spectral_density_gausswin(Ls,n,Lw_gauss);
		Tf  = sqrt(Sf);
	case {'rect'}
	%	fcut    = rectwin_cutoff_frequency(m);
%	m_gauss = fcut2Lw_gausswin(fcut);
	%^S(:,3)=normpdf(fx,0,Lw_gauss);
%	Lw = L/m
	% length of window in frequency-space
	% Ls = m*df = m/L
%	fc = rectwin_cutoff_frequency(Lw_spectral);
	Lws = m/L;
	Ls = n/L;
	% T = sqrt(S), when S is strictly positive, but this is here not the case
	[Sd,T] = sd_rectwin(Ls,n,Lws);
	end
	%fcut    = rectwin_cutoff_frequency(L/m);
	%fcut_gauss = fcut2Lw_gausswin(fcut);
	%^S(:,3)=normpdf(fx,0,Lw_gauss);
%	Sd_ = sd_rectwin(L,n,fcut);

	
	Shat = abs(fft(b)).^2;
	Rhat = ifft(Shat);
	R    = Tf.*Rhat;

	% for sanity, should not happen
	S    = real(fft(R));
	S(S<0) = 0;

	% normalize int_0^inf S df = 1
	S = 2*L*S./sum(S);
end

