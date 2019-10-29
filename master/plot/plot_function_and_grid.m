clf
subplot(3,2,1)
X=linspace(-10,10,1000);
plot(X,exp(-abs(X)),'k');
ylabel('R_1')
xlabel('\rho')
% axis equal;
hold on
%subplot(8,2,5)
 X=-10:10;
 e=1/3;
 d = X(end)/(e*(exp(e*X(end))-1));
 X = d*e*sign(X).*(exp(e*abs(X)) - 1);
% X = X(2:end-1);
 plot(X,-0.125*ones(size(X)),'k.')
ylim([-0.25 1.125])

print -deps ../../img/wf_and_grid.eps

