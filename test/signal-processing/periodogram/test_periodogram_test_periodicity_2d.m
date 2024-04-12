% Tue 20 Dec 14:34:34 CET 2022
%pp = 0:0.01:0.2;
pp = [0:0.1];
%0.05:0.05:0.95];

n = 100;
L = 1;
x = repmat(linspace(0,L,n)',1,n);
y = x';
fc = 20;

% number of bins used for smoothing
nf = 5;
% maks
bmsk = false(n,n);
m1 = n;
m2 = n;
bmsk(1:m1,1:m2) = true;
fmsk = [];
d = [];
p = [];
mu = [];
n_mc = 100;

significance_level = 0.05; 

m = 100;
p_ = zeros(m,2);
for idx=1:length(pp)
for jdx=1:m
z_periodic = sin(2*pi*fc*x);
z_random   = randn(n,n);
z_periodic = z_periodic/std(z_periodic,[],'all');
z_random   = z_random/std(z_random,[],'all');
	[idx,jdx]
	% construct a pattern with superimposed periodic frequency component
	z = pp(idx)*z_periodic + (1-pp(idx))*z_random;

	% test with mask
	[isperiodic,p_(jdx,1),stat] = periodogram_test_periodicity_2d(z, [L,L], nf, bmsk,fmsk, n_mc,significance_level);
	%p1(idx,1) = stat.p1;
	
	% test without mask (full pattern)
	[isperiodic,p_(jdx,2),stat] = periodogram_test_periodicity_2d(z, [L,L], nf, [],fmsk,[],significance_level);
	%p1(idx,2) = stat.p1;
end % for jdx
	p(idx,:) = median(p_);
	issignificant = p_<=significance_level;
	P(idx,:) = mean(issignificant);
end % for idx

clf()
subplot(2,2,1);
plot(pp,p)
xlabel('Relative magnitude of periodic component');
ylabel('Median p');
legend('Monte-Carlo','Analytic');
subplot(2,2,2);
plot(pp,P)
xlabel('Relative magnitude of periodic component');
ylabel('P(p<=a)')

