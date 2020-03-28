clear
for idx=3:16
	N(idx) = 2^idx;
	n = N(idx);
	X = pi*((0:n)/n)';
%	X = pi*((0:n)/n).^2';
%	X = pi/N(idx)*(0:N(idx))';
	f = sin(X(2:end-1));
	g = cos(X(2:end-1));
	D2_3 = fdm_d_vargrid(X, 2, 3);

	D1_5 = fdm_d_vargrid(X, 1, 5);
	D2_5 = fdm_d_vargrid(X, 2, 5);
	D3_5 = fdm_d_vargrid(X, 3, 5);
	D4_5 = fdm_d_vargrid(X, 4, 5);

%	R = 2:N(idx)-2; X = X(3:end-2);
%	R = 3:N(idx)-3; X = X(4:end-3);
	R = 4:N(idx)-4; X = X(5:end-4);

	d1f_5 = D1_5*f; d1f_5=d1f_5(R);
	d2f_5 = D2_5*f; d2f_5=d2f_5(R);
	d3f_5 = D3_5*f; d3f_5=d3f_5(R);
	d4f_5 = D4_5*f; d4f_5=d4f_5(R);
	d12f_5 = D2_5*(D1_5*f); d12f_5=d12f_5(R);
	d22f_5 = D2_5*(D2_5*f); d22f_5=d22f_5(R);

	f = f(R);
	g = g(R);

	E(idx,1) = norm(d1f_5 - g);
	E(idx,2) = norm(d2f_5 + f);
	E(idx,3) = norm(d3f_5 + g);
	E(idx,4) = norm(d4f_5 - f);
	E(idx,5) = norm(d12f_5 + g);
	E(idx,6) = norm(d22f_5 - f);

	subplot(2,3,1)
	plot(X, d1f_5 ,'.-')
	subplot(2,3,2)
	plot(X, d2f_5,'.-')
%	ylim([-1 1])
	subplot(2,3,3)
	plot(X, d3f_5,'.-')
%	ylim([-1 1])
	subplot(2,3,4)
	plot(X, d4f_5,'.-')
%	ylim([-1 1])
%pause
%	pause(1)

	%E(idx,1) = eigs(D2_3,[],1,'SM') + 1;
	%E(idx,2) = eigs(D2_5,[],1,'SM') + 1;
end
subplot(2,3,6)
loglog(N,E)

