% Thu 26 Aug 17:36:01 CEST 2021
function [c,c2] = coherence(a,b,m,flag)
	% via smoothing
	if (~flag)
	n   = floor(length(a)/m);
	%l   = length(a);
	l   = n;
	Sa  = 0;
	Sb  = 0;
	Sab = 0;
%clf
	for idx=1:m
		k   = 1+(idx-1)*n:idx*n;
		ak = a(k);
		bk = b(k);
		ak = ak -mean(ak);
		bk = bk-mean(bk);
		Fa  = fft(ak,l);
		Fb  = fft(bk,l);
		Sa  = Sa  + Fa.*conj(Fa);
		Sb  = Sb  + Fb.*conj(Fb);
		Sab = Sab + Fa.*conj(Fb);
	end
	Sa  = Sa/m;
	Sb  = Sb/m;
	Sab = Sab/m;
	else
		Fa = fft(a);
		Fb = fft(b);
		Sa  = Fa.*conj(Fa);
		Sb  = Fb.*conj(Fb);
		Sab = Fa.*conj(Fb);
		Sa  = meanfilt1(Sa,m);
		Sb  = meanfilt1(Sb,m);
		% note : smoothing must occurr before taking the absolute value
		Sab = meanfilt1((Sab),m);
	end
	c   = abs(Sab)./sqrt(Sa.*Sb);
	c2  = c.*c;
end

