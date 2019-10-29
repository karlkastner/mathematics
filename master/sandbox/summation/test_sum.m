% important : sign off summands has to change!, so no difference for norm, but
% for inner product !

%- sequential, pairwise, kahan, (sorted) 
%matlab
%intel
%atlas




% todo matrix matrix
% todo serial for comparison

function test_sum()

k=6;
t=4;

e = [0 0 0];
E_ip = zeros(k,3);
E_ip_rel = zeros(k,3);
E_n = zeros(k,3);
E_n_rel = zeros(k,3);
for idx=1:k
	n=10^idx;
	tic;
	m=0;
	while(toc < t);
		X=randn(n,1);
		Y=randn(n,1);
	
		Xs=single(X);
		Ys=single(Y);

		%S = X'*Y;
		% inner product
		S = sum_pairwise(X.*Y);
		e(1) = abs((S-Xs'*Ys));
		e(2) = abs((S-sum_pairwise(Xs.*Ys)));
		e(3) = abs((S-sum_kahan(Xs.*Ys)));
		E_ip(idx,:) = E_ip(idx,:) + e;
		E_ip_rel(idx,:) = E_ip_rel(idx,:) + e/abs(S);

		% norm ||x||^2
		S = norm(X);
		e(1) = abs(S-norm(Xs));
		e(2) = abs(S-sqrt(sum_pairwise(Xs.^2)));
		e(3) = abs(S-sqrt(sum_kahan(Xs.^2)));
		E_n(idx,:) = E_n(idx,:) + e;
		E_n_rel(idx,:) = E_n_rel(idx,:) + e/abs(S);

%		S = sum_pairwise(X);
%		es = es+abs(S-sum(Y));
%		ep = ep+abs(S-sum_pairwise(Y));
%		ek = ek+abs(S-sum_kahan(Y));

		m=m+1;
	end
	N(idx) = n;
	E_ip(idx,:)     = E_ip(idx,:)/m;
	E_ip_rel(idx,:) = E_ip_rel(idx,:)/m;
	E_n(idx,:)      = E_n(idx,:)/m;
	E_n_rel(idx,:)  = E_n_rel(idx,:)/m;
end % idx

subplot(2,2,1); loglog(N,E_ip)
title('Absolute Error of Inner-Product Calculation Depending on the Accumumaltion Algorithm');
legend('Location','NorthWest','matlab', 'pairwise', 'compensated')
grid on; set(gca,'minorgrid','none')
subplot(2,2,2); loglog(N,E_ip_rel)
grid on; set(gca,'minorgrid','none')

subplot(2,2,3); loglog(N,E_n)
title('Absolute Error of L_2-Norm Calculation Depending on the Accumumaltion Algorithm');
legend('Location','NorthWest','matlab', 'pairwise', 'compensated')
grid on; set(gca,'minorgrid','none')
subplot(2,2,4); loglog(N,E_n_rel)
grid on; set(gca,'minorgrid','none')

end % function test_sum

