% Di 12. Jan 14:27:40 CET 2016
% Fri  2 Jul 18:22:26 CEST 2021
% Karl Kastner, Berlin
%
%% estimate the spectral density nonparametrically with Bartlett's method
%
% function [S, S_std, Se, SS] = periodogram_bartlett(y,L,m,ni,pwin,subtract_mean)
% y  : data in real space
% L  : length of transect
% m  : number of segments L_seg = L/m
% ni : number of samples in output
% piwin : apply window to segments
% subtract_mean : subtract_local mean of segments
% TODO welsh sliding window
% TODO should we better normalize the individual or joint density?
% note : fourier interpolation (padding with zeros) for each segment
%        cannot be replaced by fourier interpolation of the joint density,
%        the latter differs and can have slightly negative values
function [S, S_std, Se, SS] = periodogram_bartlett(y,L,m,ni,pwin,subtract_mean)
	if (isvector(y))
		y = cvec(y);
	end

	% length of the individual segments
	n = size(y,1);
	if (m<1 || m>n)
		error('m must be larger equal 1 and smaller equal n');
	end
	np = floor(n/m);
	L_ = L*np/n;
%	L_ = L/m;
	if (nargin()<3)
		L = np;
	end
	if (nargin()<4 || isempty(ni))
		% ni = 2.^ceil(log(n)/log(2));
		ni = np;
	end
	df = 1/L;
	if (nargin>4 && ~isempty(pwin))
		w = tukeywin(ni,pwin);
		% w = w/mean(w);
	else
		w = 1;
	end
	% compute periodogram for individual segments
	siz = size(y);
	siz(1) = ni;
	S  = zeros(siz); %ni,size(y,2));
	S2 = zeros(siz); %ni,size(y,2));
%	S0 = abs(fft(y)).^2;
	if (nargout()>3)
		if (size(y,2)>1)
			error('only possible for single-vector y');
		end
		SS = zeros(ni,m);
	end
	for idx=1:m
		for jdx=1:size(y,3)
		y_(:,:,jdx) = y((idx-1)*np+(1:np),:,jdx);
		end
		if (nargin()>5 && subtract_mean)
			y_ = y_-mean(y_);
		end
 		Si = periodogram(y_,L_,ni);
		S2 = S2 + Si.*Si;
		S = S + Si;
		if (nargout>3)
			SS(:,idx) = Si;
		end 
	end
	S = S/m;
	S2 = S2/m;
	S_std = sqrt(m/(m-1).*(S2 - (S.*S)));
	Se = S_std/sqrt(m);
end

