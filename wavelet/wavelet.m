% Sun 10 Jul 13:30:12 CEST 2016
% Karl Kastner, berlin
%
%% wavelet windows
%
function [psi, win, sine] = wavelet(dt,f0,n,winstr)
	% filter length
	fs   = 1./dt;

	% avoid even length wavelet size
	%n_   = ceil(n*fs/f0);
	n_   = 2*ceil(n/2*fs/f0)+1;

	% time
	% this has indeed start at 1, not at 0 or 0.5 as otherwise the phase is 1 tap off
	t   = dt*((1:n_)'-n_/2);
	%t   = dt*((1:n_)'-1-n_/2);

	% wave : cos + 1i sin
	sine  = exp(2*pi*1i*t*f0);

	% window
	switch (lower(winstr))
	case {'rectwin'}
		win = rectwin(t);
	case {'triangular','triwin'}
		win = triwin(t);
	case {'hanning','hanwin'}
		win = hanwin(t); %,0,dt*n_/2);
	case {'kaiser', 'kaiserwin'}
		win = kaiserwin(t,0,dt*n_/2,3);
	case {'flattop','flattopwin'}
		win = flattopwin(t+dt*n_/2,0,dt*n_);
	case {'lanczos','lanczoswin'}
		win = lanczoswin(t);
	otherwise
		error('here')
	end

	% get wavelet
	psi  = win.*sine;

	% normalise
	npsi = psi'*real(sine);
	psi = (1/npsi)*psi;
end % wavelet

