%correct sign
%- implement lanczos
%- 2D - just refine one axis
% ghep with additional bc and has different first eigenvalue


%radial equation
%-1/2*(1/r^2 d/dr r^2 d/dr - l(l+1)/r^2) + V) 1/r phi = 1/r phi E

% Thu Nov  3 19:41:26 MSK 2011
% Karl Kästner, Berlin

%todo 2D, 3D
%todo, bug if L != 1

% generate convergence curve of the
function E_out = test_convergence()
opengl neverselect
tic

dimension = 2; % 3D at most 2^5
p_grid = 1;
p_fdm = 2;

fix_odd = 0;
singularity_fix = 'none'
%singularity_fix = 'clip'
%singularity_fix = 'cut_off'
%singularity_fix = 'gauss'
%singularity_fix = 'imaginary'
mode = '';
%mode = 'full';
%mode = 'GHEP';
%N=2.^(4:16)'; % 2D
%N=[2.^(4:9)'; 768]+1; % 2D
N=[2.^(4:9)']; % 2D
%N=2.^(2:4)'; % 3D
k=N(1);
c=0;
L0=10; %2.0; %3.5;
x0=0; %0.1*pi; %(1-1/sqrt(2)); %0.0;
opt = struct();
% singularity fix width
eig_mode = 'SM';
%eig_mode = 'LR';

k = 15;
kp = k-1; %6;

%Px = mu_fact;
%N = 2^16*ones(length(mu_fact),1);

%Mu = 12./N;
h = 2*L0/(N(end)+1);
%Mu = 0.1*ones(length(N),1)/h;
Mu = 12*ones(length(N),1)*2*L0/N(end);
%Mu = 1*ones(length(N),1)/sqrt(N(end));
%Mu = 12*ones(length(N),1)/N(end)^2;
Px = N;

V_cell = {};
X_cell = {};

for idx=1:length(N)
	n=N(idx);
	disp(n)

	[A X L V_] = hydrogen_boxed(n, x0, L0, dimension, p_fdm, p_grid, singularity_fix, Mu(idx));
	X_cell{idx} = X;

	if (1 == fix_odd)
		nh  = (n-1)/2;
		nnh = (size(L,1)-1)/2;
%		X(nh:nh+2)
%		pause
		X   = X([1:nh nh+2:end]);
		L   = L([1:nnh nnh+2:end], [1:nnh nnh+2:end]);
		A   = A([1:nnh nnh+2:end], [1:nnh nnh+2:end]);
		V_  = V_([1:nnh nnh+2:end], [1:nnh nnh+2:end]);
	end

	% compute eigenvalues
	switch (mode)
		case {'full'}
			E = sort(eig(full(A)));
		case { 'GHEP' }
			Vi = diag(sparse(1./diag(V_)));
			%[V E] = eigs(Vi*L+speye(size(L)),Vi,k-1,eig_mode,opt);
			%[V E] = eigs(L*Vi*Vi+Vi,Vi*Vi,k-1,eig_mode,opt);
			%[V E] = eigs(L*Vi+speye(size(L)),Vi,k-1,eig_mode,opt);
			%[V E] = eigs(Vi*L*Vi+Vi,Vi.^2,k-1,eig_mode,opt);
% del psi + 1/r psi = - lambda psi
% phi = r psi
% del 1/r phi + 1/r^2 phi = -lambda 1/r phi
				%[V E] = eigs(L*V_+ V_*V_, V_, k-1,eig_mode,opt);
			[V E] = eigs(V_*L*V_+ V_*V_*V_, V_*V_, k-1,eig_mode,opt);
			E = diag(E);
			[E ip] = sort(E);
			V = V(:,ip);
			if (1 == fix_odd)
				V = [V(1:nnh,:); zeros(1,size(V,2)); V(nnh+1:end,:)];
			end
			[X V] = attach_boundary_value(X_cell{idx},V,L0,dimension);
			X_cell{idx} = X;
			E_cell{idx} = E;
			V_cell{idx} = V;
		otherwise
			[V E] = eigs(A,k-1,eig_mode,opt);
			E = diag(E);
			[E ip] = sort(E);
			V = V(:,ip);
			if (1 == fix_odd)
				V = [V(1:nnh,:); zeros(1,size(V,2)); V(nnh+1:end,:)];
			end
			[X V] = attach_boundary_value(X_cell{idx},V,L0,dimension);
			E_cell{idx} = E;
			V_cell{idx} = V;
			X_cell{idx} = X;
	end % switch
	
	En(:,idx) = E(end-k+2:end);
end % for idx

M = V_*L*V_+ V_*V_*V_; sum(sum((M'-M).^2))
%V_ = speye(size(V,1));V_(:,1:size(V,2)) = V;cond(V), pause

% convert eigenvalues into energy levels
Eh=1;%Eh = 4.35974394e-18; % hartree in J
En = Eh*En;

f=1; hf = figure(f);% set(h,'Units','normalized','Position',[0 0 1 1])
for idx=1:kp %length(En(:,1))
	%subplot(4,4,idx);
	subplot(ceil(sqrt(kp)),ceil(kp/ceil(sqrt(kp))),idx);
	semilogx(Mu, En(idx,:));
end
f=f+1; figure(f)

% relative error, TODO, why not absolute
na = norm(En(:,end));
errA = abs(En(:,1:end-1) - En(:,end)*ones(1,size(En,2)-1))/na;

% norm of the error
normerrA = sqrt(sum(errA.^2))';

disp('eigenvalues')
disp(En(:,end)')
E_out = En(:,end)';
disp('error')
disp(errA(:,end)')

% rate of convergence
C=log([errA(:,1)]./[errA(:,end)]) / (log(N(end-1)/N(1)));

% find good and bad eigenvalue, e.g. those who converged at least with rate 1 and those who did not
good = find(C(:,1) > 1);
bad = find(C(:,1) <= 1 );

normerrA_good = sqrt(sum(errA(good,:).^2,1))';
normerrA_bad = sqrt(sum(errA(bad,:).^2,1))';

f=f+1; figure(f)
subplot(3,2,1);
loglog(Px(1:end-1),[normerrA normerrA_good normerrA_bad]);
grid on; set(gca,'minorgrid','none');

subplot(3,2,3);
plot(En(:,2:end)','.');

f=f+1; figure(f);
subplot(2,2,1);
loglog(Px(1:end-1),errA,'.-');
grid on; set(gca,'minorgrid','none');
title(['Convergence of individual eigenvalues L_0=' num2str(L0) ' x_0=' num2str(x0) ' fix=' singularity_fix])
ylim([1e-12 1]);
xlabel('number of grid points N');
ylabel('estimated error |\lambda_{N} - \lambda|')
s = sprintf('convergence-invdividual-l0-%1.2f-x0-%1.2f-fix-%s.eps',L0,x0,singularity_fix)
print('-deps',[ '../img/' s])

% probability density function
%for jdx=1:length(V_cell)
for jdx=length(V_cell)-1:length(V_cell)
E = E_cell{jdx};
V = V_cell{jdx};
X = X_cell{jdx};
f=f+1; figure(f)
clf
W = real(V.^2);
for idx=1:kp
	%subplot(kp,2,2*idx-1);
	subplot(ceil(sqrt(kp)),ceil(kp/ceil(sqrt(kp))),idx);
	plot_wavefunction(X,W(:,end-idx+1),dimension)
	title(E(end-idx+1,end))
end % for idx
end

% spatial distribution of the error

for jdx=length(V_cell)-1:length(V_cell)-1
X_half = X_cell{jdx};
V_half = V_cell{jdx};
f=f+1; figure(f)
clf
for idx=1:kp
	% interpolate
	switch dimension
		case {1}
			V_half_(:,idx) = interp1(X_half,V_half(:,idx),X);
		case {2}
			n = sqrt(size(V_half(:,1),1));
			W = reshape(V_half(:,idx),n,n);
			[XX YY] = meshgrid(X,X);
			[XX_half YY_half] = meshgrid(X_half,X_half);
			W = interp2(XX_half,YY_half,W,XX,YY);
			n = length(X);
			W = reshape(W,n^2,1);
			V_half_(:,idx) = W;
	end
	% renorm
	V_half_(:,idx) = V_half_(:,idx)/norm(V_half_(:,idx));
end

Werr = V.^2 - V_half_.^2;
for idx=1:kp
	%subplot(kp,2,2*idx-1);
	subplot(ceil(sqrt(kp)),ceil(kp/ceil(sqrt(kp))),idx);
	plot_wavefunction(X,Werr(:,end-idx+1),dimension)
	title(E(end-idx+1,end))
end % for idx
end

toc
end % function test_convergence

