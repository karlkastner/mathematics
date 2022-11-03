y = randn(1e4,100);
L = 1;
%ni =

se = [NaN,NaN];
for idx = 2:10
	m = idx;
	[S, S_std, Se] = periodogram_bartlett(y,L,m); %,ni,pwin,subtract_mean)
	se(idx,1) = rms(flat(Se));
	Se = std(S,[],2);
	se(idx,2) = rms(Se);
end
plot(se)

