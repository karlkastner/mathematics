% 2015-11-09 12:23:45.782162063 +0100
% Karl Kastner, Berlin
%
%% continuous fourier transformation of y
%% (not discrete fourier transformation dft/fft)
%%
%% input:
%%	b : data sampled at equal intervals
%%	T : length of data in time or space, i.e. position of last sample if
%%	    position of first sample is 0
%%       T_max : maximum period to include
%% 
%% output :
%%	A  : fourier matrix
%%       p  : fourier transformation of b
%%	tt : TODO 
%
% TODO allow user to choose arbitrary base period
% TODO split into matrix setup and transformation
function [A, p, tt] = fourier_transform(y,T,T_max)
	% create the Fourier matrix
	n = length(y);
	t = T/n*(0:n-1)';
	m = floor(n/2);
	% only include durations shorter than T_max
%	m_min = round(T/T_max);
	m_ = round(T_max*n/T)
	S = zeros(n,m_);
	C = zeros(n,m_);
	for idx=1:m_; %m_min:m
		dt = T/n*idx;
		f = 1/dt;
		S(:,idx) = sin(2*pi*f*t);
		C(:,idx) = cos(2*pi*f*t);
		%S(:,idx-m_min+1) = sin(2*pi*t*idx/T);
		%C(:,idx-m_min+1) = cos(2*pi*t*idx/T);
		%S(:,idx-m_min+1) = sin(2*pi*t*idx/T);
		%C(:,idx-m_min+1) = cos(2*pi*t*idx/T);
	end
%	A  = [ones(n,1) S C];
	A  = [S, C];
	tt = T./[1:m, 2:m];
	% regress the coefficients
	p = A \ y;
end % fourier analyis

