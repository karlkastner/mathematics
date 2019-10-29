% 2012 Jan  2  2012
% Karl Kästner, Berlin

% http://panda.unm.edu/Courses/Finley/P262/Hydrogen/WaveFcns.html

n = 1000;
L0 = 10;
X = L0*((0:n)/(n+1) - 0.5);
r = abs(X);
a0 = 1;
%Phi = 1/sqrt(pi)*(1/a0)^(3/2)*exp(-R/(2*a0));
% radial wavefunction
R = 2*(1/a0)^(3/2)*exp(-r/a0);
D = (r.*R).^2;
D = D/norm(D)
figure(1)
subplot(4,2,1)
plot(X,D,'k','Linewidth',2)
xlabel('x in a_0')
title('Radial Probabilty Distribution r^2 R^2 of the Hydrogen Ground State')
axis tight
print -deps ../img/hydrogen-radial-probability.eps

figure(2)
clf
subplot(4,2,1)
%R200 = (1/(2*a0))^(3/2)*(2 - abs(X)/a0).*exp(-abs(X)/(2*a0));
R200 = [];
plot(X, [R' R200'], 'k', 'Linewidth',2)
xlabel('r [a_0]')
ylabel('\Psi(r)')
title('Radial Wave Function of the Hydrogen Ground State')
axis tight
print -deps ../img/hydrogen-radial-wf.eps
system('epstopdf ../img/hydrogen-radial-wf.eps');

figure(3)
clf
for n=1:2
L=1;
N=100;
X = (0:N)'*L/N;
k=n*pi/L;
psi = sqrt(2/L)*sin(k*X);
subplot(4,2,1); hold on
plot(X,psi,'k', 'Linewidth',2);
title('Wave Function of the Lowest two States of a Particle in a Box')
xlabel('x')
ylabel('\Psi')
end
print -deps ../img/box-wf.eps

figure(4)
clf
subplot(2,2,1)
n=1:16;
E=-13.6*1./n;
for idx=1:length(n)
	%plot(n,E,'.')
	line([n(idx) n(idx)], [E(idx) 0],'color',[0 0 0],'linewidth',2); hold on
end
title('Energy Levels of the Unconfined Hydrogen Atom')
xlabel('principal quantum number n')
ylabel('E in eV')
xlim([0 n(end)+1])
ylim([-15 0])
grid on
print -deps ../img/hydrogen-energy-levels.eps
system('epstopdf ../img/hydrogen-energy-levels.eps')

figure(5)
clf
subplot(2,2,1)
n=1:16;
E=n.^2;
for idx=1:length(n)
	%plot(n,E,'.')
	line([n(idx) n(idx)], [0 E(idx)],'color',[0 0 0],'linewidth',2); hold on
end
title('Energy Levels of a single Particle in the 1D unit Box')
xlabel('energy level k')
ylabel('E in eV')
xlim([0 n(end)+1])
%ylim([-15 1])
grid on
print -deps ../img/box-energy-levels.eps

