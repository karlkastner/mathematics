% 2013-08-11 17:44:52.000000000 +0200
% Karl Kastner, Berlin

% scaling : amplitude is attenuated
% cross talk / non orthogonality : side lobes
% phase ?

n = 1e3;
m = 100;

T = [1 0.5 0.25];
C_ = [];

%wf = 'haar';
wf = 'morl';

% morl resolves frequencies better and has almost no side lobes

for idx=1:length(T)
	n=1e3;
	t = linspace(0,10,n);
	f = logspace(0,log(n)/log(10),m);
	A = sin(2*pi*t/T(idx)) + sin(2*2*pi*t/T(idx));
%	A = sin(2*pi*t/T(idx));
	C = cwt(A,f,wf);
	C_ = [C_; sum(abs(C'))./sqrt(f)];	% 1/sqrt(f) scaling only correct for haart wavelet
end
	%imagesc(abs(C));
	semilogx(f,C_')
	colormap gray;

