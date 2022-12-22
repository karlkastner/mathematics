% Fri 23 Sep 16:22:32 CEST 2022
% fcut : frequency, at which the window reaches 0.5
function w = gausswin(f,fcut,df)
	s = fcut/sqrt(-log(0.25));
	w = normpdf(f,0,s)/normpdf(0,0,s);
	w = w./(sum(w).*df);
end

