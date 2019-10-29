function test_regularisation
c = {'r','g','b','k'}
N = 2.^(8:11);
clf
V_ = [];
for idx=1:length(N)
n=N(idx);
x0 = 0;
x = (-n:n)'/n;
x = (-n+0.5:n-0.5)'/n;
Vs = potential(x);
V = 1./abs(x(2:end-1));
%plot(x(2:end-1),[V Vs]); hold on
plot(x(2:end-1),[V-Vs],c{idx}); hold on
%V_(:,idx) = eigs(diag(sparse(Vs)),6,'LM')
V_(:,idx) = eigs(diag(sparse(V)),6,'LM');
end
V_
end

function V = potential(x)
	% todo - this is only second order accurate
	xs = sqrt(abs(x-x0));
	h = 1/(n;
	Vs = (( xs(3:end)-xs(1:end-2) )/h).^2;
end

