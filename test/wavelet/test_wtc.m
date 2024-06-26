% Thu Dec 26 16:38:43 WIB 2013
% Karl Kastner, Berlin
%
% test of the wtc package by Grinsted
function test_wtc()

addpath([ROOTFOLDER,'src/wtc']);

dt = 0.01;
t = 1:dt:10;
f = 1;

D = sin(2*pi*f*t);
%wt(D);
[wave period scale coi] = wt(D);
wave_ = cwt(D,period,'morl');
period = period*dt;
subplot(1,2,1)
imagesc(t,period,abs(wave))
subplot(1,2,2)
imagesc(t,period,abs(wave_))

pause
%plot(coi)
%dt*period
%coi
%pause
f_ = 2.^(-10:10);
F = f*logspace(-2,1,51);

for jdx=1:length(f_) 

D = sin(2*pi*f_(jdx)*t);
for idx=1:length(F)
%	T_min = dt/F(idx);
%	T_max = 
	% scales are actually periods not frequencies
	[wave period scale coi] = wt(D,'s0',F(idx)/dt, 'MaxScale',1.1*F(idx)/dt);
	period_ = period;
	period = period(1);
	fdx = find(coi > period);
	l = length(fdx);
	wave = wave(1,:)/sqrt(scale(1));
	A(idx,jdx) = norm(wave(fdx))/sqrt(l);
	
	
end % for idx

end
%semilogy(F,A);
semilogx(F,A);
grid on

% what do period and scale mean?

end

