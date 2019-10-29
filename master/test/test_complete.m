%1D energy levels
%-> raw values | eV
%-> convergence with n
%-> influence of cut-off  [0.1 1 10 100] h_min
%-> influence of boundary [1 10 100]
%-> influence of location inside the domain

function test_complete()

%
% 1D
%

% 1D Convergence
if (0)
clear;
n = ceil(logspace(log(2^4)/log(10),log(2^16)/log(10),10));
dimension = 1;
singularity_fix = 'cut_off';
L0 = 2;
x0 = 0;
k = 16;
p_fdm = 4;
p_grid = 1;
a = [0.1 1 10 100]/(n(end)+1);

for adx=1:length(a)
	adx
 for ndx=1:length(n)
	ndx
        A= hydrogen_boxed(n(ndx), x0, L0, dimension, p_fdm, p_grid, singularity_fix, a(adx));
        E(:,ndx) = sort(eigs(A,k,'SM'));
 end % ndx
 Err=E(:,1:end-1)-E(:,end)*ones(1,length(n)-1);
 nErr(adx,:) = sqrt(sum(Err.*Err));
end % adx
subplot(2,2,1);
loglog(n(1:end-1),nErr,'.-');
title('1D Eigenvalues Convergence'); % todo - low mu, medium mu, high mu
xlabel('n');
ylabel('Eigenvalue');
grid on;
print -depsc ../img/1D-convergence.eps
end


% 1D Cut-Off influence
if (0)
clear;
n = 2^16;
dimension = 1;
singularity_fix = 'cut_off';
L0 = 2;
x0 = 0;
k = 16;
p_fdm = 4;
p_grid = 1;
a = logspace(1,2,10)/(n+1);

for adx=1:length(a)
	adx
        A= hydrogen_boxed(n, x0, L0, dimension, p_fdm, p_grid, singularity_fix, a(adx));
        E(:,adx) = sort(eigs(A,k,'SM'));
end
subplot(2,2,1);
semilogx(a*(n+1),E,'.')
title('1D Eigenvalues Depending on Cut-Off Potential V=1/|x+a|')
xlabel('a/h')
ylabel('Eigenvalue')
grid on
print -depsc ../img/1D-cut-off-influence.eps
end

% 1D influence of box-size
if (0)
clear;
n = 2^16;
dimension = 1;
singularity_fix = 'cut_off';
L0 = logspace(0,2,10);
x0 = 0;
k = 16;
p_fdm = 4;
p_grid = 1;
a = 10/(n+1);

for ldx=1:length(L0)
	ldx
        A= hydrogen_boxed(n, x0, L0(ldx), dimension, p_fdm, p_grid, singularity_fix, a);
        E(:,ldx) = sort(eigs(A,k,'SM'));
end
subplot(2,2,1);
semilogx(L0, E, '.')
title('1D Eigenvalues Levels Depending on Box-Size')
xlabel('L0')
ylabel('Eigenvalue')
grid on
print -depsc ../img/1D-box-size-influence.eps
end

% 1D influence of proton-location
if (0)
clear;
n = 2^16;
dimension = 1;
singularity_fix = 'cut_off';
L0 = 2;
x0 = linspace(0,L0,10);
k = 16;
p_fdm = 4;
p_grid = 1;
a = 10/(n+1);

for ldx=1:length(x0)
	ldx
	x0(ldx)
        A= hydrogen_boxed(n, x0(ldx), L0, dimension, p_fdm, p_grid, singularity_fix, a);
        E(:,ldx) = sort(eigs(A,k,'SM'));
end
subplot(2,2,1);
plot(x0/L0, E, '.')
title('1D Eigenvalues Levels Depending on Proton-location')
xlabel('x0/L0')
ylabel('Eigenvalue')
grid on
print -depsc ../img/1D-proton-location-influence.eps
end

%
% 2D
%

% 2D Convergence
if (0)
clear;
n = ceil(logspace(log(2^4)/log(10),log(2^8)/log(10),10));
dimension = 2;
singularity_fix = 'cut_off';
L0 = 2;
x0 = 0;
k = 16;
p_fdm = 4;
p_grid = 1;
a = [1 10 100]/(n(end)+1);

for adx=1:length(a)
	adx
 for ndx=1:length(n)
	n(ndx)
        A= hydrogen_boxed(n(ndx), x0, L0, dimension, p_fdm, p_grid, singularity_fix, a(adx));
        E(:,ndx) = sort(eigs(A,k,'SM'));
 end % ndx
 Err=E(:,1:end-1)-E(:,end)*ones(1,length(n)-1);
 nErr(adx,:) = sqrt(sum(Err.*Err));
end % adx
subplot(2,2,1);
loglog(n(1:end-1),nErr,'.-');
title('2D Eigenvalues Convergence'); % todo - low mu, medium mu, high mu
xlabel('n');
ylabel('Eigenvalue');
grid on
print -depsc ../img/2D-convergence.eps
end

% 2D Cut-Off influence
if (0)
clear;
n = 2^8;
dimension = 2;
singularity_fix = 'cut_off';
L0 = 2;
x0 = 0;
k = 16;
p_fdm = 4;
p_grid = 1;
a = logspace(-1,2,10)/(n+1);

for adx=1:length(a)
	adx
        A= hydrogen_boxed(n, x0, L0, dimension, p_fdm, p_grid, singularity_fix, a(adx));
        E(:,adx) = sort(eigs(A,k,'SM'));
end
subplot(2,2,1);
semilogx(a*(n+1),E,'.')
title('2D Eigenvalues Depending on Cut-Off Potential V=1/|x+a|')
xlabel('a/h')
ylabel('Eigenvalue')
grid on
print -depsc ../img/2D-cut-off-influence.eps
end

% 2D influence of box-size
if (0)
clear;
n = 2^8;
dimension = 2;
singularity_fix = 'cut_off';
L0 = logspace(0,2,10);
x0 = 0;
k = 16;
p_fdm = 4;
p_grid = 1;
a = 10/(n+1);

for ldx=1:length(L0)
	ldx
        A= hydrogen_boxed(n, x0, L0(ldx), dimension, p_fdm, p_grid, singularity_fix, a);
        E(:,ldx) = sort(eigs(A,k,'SM'));
end
subplot(2,2,1);
semilogx(L0, E, '.')
title('2D Eigenvalues Levels Depending on Box-Size')
xlabel('L0')
ylabel('Eigenvalue')
grid on
print -depsc ../img/2D-influence-of-box-size.eps
end

% 2D influence of proton-location
if (1)
clear;
n = 2^8;
dimension = 2;
singularity_fix = 'cut_off';
L0 = 2;
x0 = linspace(0,L0,10);
k = 16;
p_fdm = 4;
p_grid = 1;
a = 10/(n+1);

for ldx=1:length(x0)
	ldx
	x0(ldx)
        A= hydrogen_boxed(n, x0(ldx), L0, dimension, p_fdm, p_grid, singularity_fix, a);
        E(:,ldx) = sort(eigs(A,k,'SM'));
end
subplot(2,2,1);
plot(x0/L0, E, '.');
title('2D Eigenvalues Levels Depending on Proton-location');
xlabel('x0/L0 axis parallel');
ylabel('Eigenvalue');
grid on;
print -depsc ../img/2D-influence-of-proton-location.eps
end

%
% 3D
%

if (1)
clear;
n = unique(ceil(logspace(log(5)/log(10),log(2^5)/log(10),10)));
dimension = 3;
singularity_fix = 'cut_off';
L0 = 2;
x0 = 0;
k = 16;
p_fdm = 4;
p_grid = 1;
a = [1 10 100]/(n(end)+1);

for adx=1:length(a)
	adx
 for ndx=1:length(n)
	ndx
	n(ndx)
        A= hydrogen_boxed(n(ndx), x0, L0, dimension, p_fdm, p_grid, singularity_fix, a(adx));
        E(:,ndx) = sort(eigs(A,k,'SM'));
 end % ndx
 Err=E(:,1:end-1)-E(:,end)*ones(1,length(n)-1);
 nErr(adx,:) = sqrt(sum(Err.*Err));
end % adx
subplot(2,2,1);
loglog(n(1:end-1),nErr,'.-')
title('2D Eigenvalues Convergence'); % todo - low mu, medium mu, high mu
xlabel('n')
ylabel('Eigenvalue')
grid on
print -depsc ../img/3D-convergence.eps
end

% 3D Cut-Off influence
if (1)
clear;
n = 2^5;
dimension = 3;
singularity_fix = 'cut_off';
L0 = 2;
x0 = 0;
k = 16;
p_fdm = 4;
p_grid = 1;
a = logspace(-1,2,10)/(n+1);

for adx=1:length(a)
	adx
        A= hydrogen_boxed(n, x0, L0, dimension, p_fdm, p_grid, singularity_fix, a(adx));
        E(:,adx) = sort(eigs(A,k,'SM'));
end
subplot(2,2,1);
semilogx(a*(n+1),E,'.');
title('3D Eigenvalues Depending on Cut-Off Potential V=1/|x+a|');
xlabel('a/h')';
ylabel('Eigenvalue');
grid on;
print -depsc ../img/3D-influence-of-cut-off.eps
end

% 3D influence of box-size
if (1)
clear;
n = 2^5;
dimension = 3;
singularity_fix = 'cut_off';
L0 = logspace(0,2,10);
x0 = 0;
k = 16;
p_fdm = 4;
p_grid = 1;
a = 10/(n+1);

for ldx=1:length(L0)
	ldx
        A= hydrogen_boxed(n, x0, L0(ldx), dimension, p_fdm, p_grid, singularity_fix, a);
        E(:,ldx) = sort(eigs(A,k,'SM'));
end
subplot(2,2,1);
semilogx(L0, E, '.')
title('3D Eigenvalues Levels Depending on Box-Size')
xlabel('L0')
ylabel('Eigenvalue')
grid on
print -depsc ../img/3D-influence-of-box-size.eps
end

% 3D influence of proton-location
if (1)
clear;
n = 2^5;
dimension = 3;
singularity_fix = 'cut_off';
L0 = 2;
x0 = linspace(0,L0,10);
k = 16;
p_fdm = 4;
p_grid = 1;
a = 10/(n+1);

for ldx=1:length(x0)
	ldx
	x0(ldx)
        A= hydrogen_boxed(n, x0(ldx), L0, dimension, p_fdm, p_grid, singularity_fix, a);
        E(:,ldx) = sort(eigs(A,k,'SM'));
end
subplot(2,2,1);
plot(x0/L0, E, '.');
title('3D Eigenvalues Levels Depending on Proton-location');
xlabel('x0/L0 axis parallel');
ylabel('Eigenvalue');
grid on;
print -depsc ../img/3D-influence-of-proton-location.eps
end

end

