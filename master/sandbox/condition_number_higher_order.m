% Mon Feb 11 00:12:25 MSK 2013
% Karl Kästner, Berlin

n=1000;
k=n;
c = [1 -1/12 1/90 -1/560 1/3150];
A0 = (n+1)^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
A=0.*speye(n);
for idx=1:length(c);
	A = A + c(idx)*A0^(idx)*1/(n+1)^(2*(idx-1));
	max(sum(abs(A)))/(n+1)^2
	E(:,idx)=eigs(A,[],k,'LM')
end
plot(-flipud(E)/(n+1)^2)
E(1,:)./(n+1)^2

