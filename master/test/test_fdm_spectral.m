% Sat Mar 24 04:15:48 MSK 2012
% Karl Kästner, Berlin

% TODO rename - this is not a test case !!!


clf
h=2.^(-7:0);
k=logspace(0,3,1e5);
subplot(2,2,1)
df = sin(k'*h) ./ (ones(length(k),1)*h) .* (k'*h < pi);
loglog(k, [k' df],'k');
for idx=1:length(h)
	%text(2.^(1:8), ones(size(h))', sprintf('2^{%d}\n',log(h')/log(2)))
	text(0.55*pi*2.^(8-idx), 2e-1, sprintf('h=2^{%d}\n',log(h(idx)')/log(2)),'fontsize',6)
end
%legend('analytic',num2str(log(h')/log(2)))
xlabel('k');
ylabel('$\frac{1}{f(k)} f_h''(k)$','interpreter','latex')
%title('First Derivative of f(x) = e^{ikx}');
ylim([1e-1 5e2]);
grid on;
set(gca,'minorgrid','none')

subplot(2,2,2)
ddf = 2*( cos(k'*h) - 1) ./ (ones(length(k),1)*h.^2) .* (k'*h < 2*pi );
loglog(k, -[-k'.^2 ddf],'k');
for idx=1:length(h)
	%text(2.^(1:8), ones(size(h))', sprintf('2^{%d}\n',log(h')/log(2)))
	text(0.55*2*pi*2.^(8-idx), 3e-2, sprintf('h=2^{%d}\n',log(h(idx)')/log(2)),'fontsize',6)
end
xlabel('k');
ylabel('$\frac{1}{f(k)} f_h''''(k)$','interpreter','latex')
ylim([1e-2 1e5]);
grid on;
set(gca,'minorgrid','none')
%plot(k,ddf)

print -deps ../../img/fdm_spectral_error.eps

