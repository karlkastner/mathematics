% Thu Dec 26 12:32:51 WIB 2013
% Karl Kastner, Berlin
%
%% continuous wavelet transform 
%% follows "The Illustrated Wavelet Transform Handbook: Introductory Theory and ..."
function [A, P] = continuous_wavelet_transforms(t,y,f)
	% TODO check for homogeneously spaced intervals
	% (take dt as parameter)
	dt = t(2)-t(1);
	% create windows
	% TODO no magic numbers
	c1 = 3;
	c2 = 1/c1;
	t_ = (-c1/f):dt:(c1/f);
	%t_ = (-c1/f):dt:(c1/f);
	s = 2*pi*f*1;
	kappa = 0;
%	w = sqrt(1 + exp(-s*s) - 2*exp(-0.75*s*s)) ...
%            * pi^-0.25*exp(-0.5*(2*pi*f*t_).^2).*(exp(1i*s*(2*pi*f)*t_) - kappa);
%	w = sqrt(1 + exp(-s*s) - 2*exp(-0.75*s*s)) ...
%           * pi^-0.25*exp(-0.5*(2*pi*f*t_).^2).*(exp(1i*(2*pi*f)*t_) - kappa);
%	w = sqrt(2*pi*s/dt)*pi^-0.25*exp(-(s*2*pi*f*t_ - 2*pi*f).^2);
%	w = pi.^-0.25*exp(1i*2*pi*f*t_ - 0.5*t_.*t_);
%	w = w / pi^2;
%	% standard morlet wavelet
%	omega0 = 5;
%	w = pi^-0.25 exp(-0.5*t_.^2 + 1i*omega0*t_);
%	w = exp(-0.5*s*t_.^2).*(exp(1i*2*pi*f*t_) - exp(-0.5*(2*pi*f)^2/s));
%	w = w/norm(w);
	% solve for central frequency
	omega = 2*pi*f;
	sigma = fzero( @(sigma) (omega - sigma).^2 - 1 - (omega^2 - 1)*exp(-sigma*omega), 1 );
	% construct the wavelet
	f_ = f;
	w = exp(-0.5*(f_*t_).^2).*(exp(1i*sigma*t_));
	w = w/norm(w);


	ws = real(w);
	wc = imag(w);
%	plot(t_, [ws' wc']);
%	f
%	title(1/f)
%	 pause;

%	w  = normpdf(t_,0,c2*2*pi/f);
%	ws = sin(2*pi*f*t_).*w;
%	wc = cos(2*pi*f*t_).*w;
%	% normalise
%	ws = ws / norm(ws);
%	wc = wc / norm(wc);

	% TODO pad nans
	% convolve
	cs = conv(y, ws, 'same');
	cc = conv(y, wc, 'same');
	% extract amplitude
	A  = sqrt(cs.*cs + cc.*cc);
	% extract phase
	P  = atan2(cs,cc);
end % wa

