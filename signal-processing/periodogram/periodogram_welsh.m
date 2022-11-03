% Fri 29 Apr 18:39:44 CEST 2022
function S = periodogram_welsh(x,L,m)
	S = 0;
	n = length(x);
	% double overlap
	p = 2;
	ni  = round(n/m);
	nie = round(p*n/m);
	% periodic extension
	x = [x;x];
	for idx=1:m
		id=(idx-1)*ni+(1:nie);
		x_ = x(id);
		S= S+abs(fft(x(id),n)).^2;
	end
	S = S/m;
	df = 1/L;
	S = 2*S/(sum(S)*df);
	S = cvec(S);
end

