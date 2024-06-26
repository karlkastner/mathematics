% Thu Dec 26 12:42:08 WIB 2013
% Karl Kastner, Berlin

function test_wavelet_analysis()

D=load([ROOTFOLDER,'src/tide/table-pontianak/level_.out']);
D = D(:,end-1);
t = (1:length(D))'/48;
f = 1; % diurnal tide
f = 2;

[A P] = wavelet_analysis(t, D, f);

% pure sine frequency response
D = sin(2*pi*f*t);
F = logspace(-2,2,15);
for idx=1:length(F)
	[A P] = wavelet_analysis(t, D, f*F(idx));
	R(idx) = max(A);
end
figure(1)
loglog(F,R);

figure(2)
% pure sine phase response
P_ = linspace(0,pi,12);
for idx=1:length(P_)
	subplot(3,4,idx)
	D = sin(2*pi*f*t + P_(idx));
	[A P] = wavelet_analysis(t, D, f);
	plot(P,'-')
	xlim([0 1000])
	title(P_(idx))
	%pause(1)
end

%plot([D A],'.')

end

