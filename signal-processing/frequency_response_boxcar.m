% 2015-07-19 19:37:41.746503824 +0200
%% frquency response of a boxcar filter
function H = frequency_response_boxcar(L,f)
	H = bsxfun(@times,rvec(1./L),cvec((1-exp(-1i*2*pi*f*L))./(1-exp(-1i*2*pi*f))));
end

