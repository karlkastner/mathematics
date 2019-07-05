% Thu Jul 25 20:47:58 UTC 2013
% Karl Kästner, Berlin
%
%% continuous fourier transform
%% as of time of implmentation, the matlab interal cwt is affected by
%% serious round-off errors and has issues with the scaling,
%% which is not the case here
%
% matlab internal function misses the sigma parameter
% matlab internal function exposes round off error
% TODO use exact transform and inverse (scalefactor!!!)
% low frequencies are underestimated, if data series is not long enough, especially at the rim
% TODO cone of influence (esp when using the max)
function c = cwt_man(x,scale,sigma)
	c = zeros(length(scale),length(x));
	for idx=1:length(scale)
		% generate a wavelet
		s = scale(idx);
		s = scale(end)/scale(idx);
% TODO, t can be short cut to -5 sigma s 
		t = ((1:length(x))-length(x)/2)/length(x);
%		length(t)
%		fdx = find(t > -5/s & t < 5/s);
%		[max(t) 5/s]
%		t = t(fdx);
%		length(t)		
%		t =  max(-length(x)/2,-5*sigma*s):min(5*sigma*s,length(x)/2);
%		[length(t) length(t_)]
		w = exp(-0.5*s^2*t.^2).*cos(sigma*s*t);
		C = 1/norm(w);
%		C = %1/norm(w.*t);
%		C = sigma^2*(s)./( sum( w.^2 .* 1./abs(t) ) );

		w = C*w;
%		C = 2*pi*sum(w./abs(t))*(t(2)-t(1))
%		C = 1/s*scale(end)*sqrt(2);
%		C = 1;
		c(idx,:) = conv(x,w,'same');
	end
end % cwt_man()

