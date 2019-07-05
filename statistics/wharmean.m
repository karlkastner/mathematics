% 2017-02-16 10:57:57.150396861 +0100
%% weighted harmonic mean
function mu = wharmean(w,x)
	w = w/sum(w);
	mu = 1./sum(w./x);
end

