
figure(1)
clf
subplot(1.5,1.5,1)
path(path,'../')
path(path,'../dat')
%eval('fdm_1d_dat')
%nErr = relerr([E([2 3],1:end) [-0.5; -0.5]]);
%loglog(N,nErr(1:end-1),'-', 'Linewidth',2,'color',[1   0 0]); hold on
%eval('fem_1d_dat')
%nErr = relerr([E([2 3],1:end) [-0.5; -0.5]]);
%loglog(N,nErr(1:end-1),'--', 'Linewidth',2,'color',[1   0 0]); hold on

	%loglog(N,relerr(E(2:2:end,:)),'r-', 'Linewidth',2); hold on
	%loglog(N,relerr(E(3:2:end,:)),'r--', 'Linewidth',2); hold on
eval('fdm_2d_dat')
nErr = relerr([E([1 2],1:size(E,2)) [-2; -2/3^2]]);
loglog(N,nErr(1:end-1),'-', 'Linewidth',2,'color',[0   0.75 0]); hold on
eval('fem_2d_dat')
nErr = relerr([E([1 2],1:end) [-2; -2/3^2]]);
loglog(N,nErr(1:end-1),'--', 'Linewidth',2,'color',[  0 0.75 0]); hold on
	%loglog(N,relerr(E([1 4 9],:)),'g--', 'Linewidth',2); hold on
	%loglog(N,relerr(E([2 3 5 6 7 8 10],:)),'g-', 'Linewidth',2); hold on
eval('fdm_3d_dat')
nErr = relerr([E([1 2],1:end) [-0.5; -0.5/4]]);
loglog(N,nErr(1:end-1),'-', 'Linewidth',2, 'color', [0 0 1]); hold on
eval('fem_3d_dat')
nErr = relerr([E([1 2],1:end) [-0.5; -0.5/4]]);
loglog(N,nErr(1:end-1),'--', 'Linewidth',2, 'color', [0 0 1]); hold on
	%loglog(N,relerr(E),'b-', 'Linewidth',2); hold on
	%loglog(N,sqrt((E(2,:)+0.5/4).^2)*4/0.5,'b'); hold on


eval('harm_fdm_1d')
nErr = relerr([E([1 2],1:end) [1; 3]]);
loglog(N,nErr(1:end-1),'-', 'Linewidth',2,'color',[0.5   0.5 0.5]); hold on

axis([3e0 1e3 1e-2 3e0])
%legend('Location','SouthWest','1D even','1D odd', '2D \lambda 1,4,9', '2D others', '3D all')
legend('Location','SouthWest','FDM 1D', 'FEM 1D', 'FDM 2D', 'FEM 2D', 'FDM 3D', 'FEM 3D');
grid on; set(gca,'minorgrid','none')

ylabel('$$\Big | \frac{\lambda_n \mbox{-}  \lambda}{\lambda} \Big |$$','interpreter','latex','Rotation',0)
xlabel('number of grid points per axis')

print -depsc ../../img/h_convergence.eps
title('Relative Discretisation Error - Curse of Dicontinuous Derivatives')

figure(2);
clf;
subplot(1.25,1.25,1);
% {
%eval('fdm_1d_vargrid');
%nErr = relerr([E([2 3],1:end) [-0.5; -0.5]]);
%loglog(N,nErr(1:end-1),'-', 'Linewidth',2,'color',[1   0 0]); hold on
 %loglog(N,abs(E(2,:)+0.5),'r-', 'Linewidth',2,'color',[0.75 0 0]); hold on
 %loglog(N,abs(E(3,:)+0.5),'--', 'Linewidth',2,'color',[0.75 0 0]); hold on
%eval('fdm_1d_dat');
%nErr = relerr([E([2 3],1:end) [-0.5; -0.5]]);
%loglog(N,nErr(1:end-1),'--', 'Linewidth',2,'color',[1   0.5 0.5]); hold on
 %loglog(N,abs(E(2,:)+0.5),'r-', 'Linewidth',2,'color',[1 0.25 0.25]); hold on
 %loglog(N,abs(E(3,:)+0.5),'--', 'Linewidth',2,'color',[1 0.25 0.25]); hold on
% }

eval('fdm_2d_vargrid_dat_new');
nErr = relerr([E([1 2],2:end) [-2; -2/3^2]]);
N=N(2:end);
loglog(N,nErr(1:end-1),'-', 'Linewidth',2,'color',[0   0.75 0]); hold on
eval('fdm_2d_dat');
nErr = relerr([E([1 2],:) [-2; -2/3^2]]);
loglog(N,nErr(1:end-1),'--', 'Linewidth',2,'color',[0.375   0.75 0.375]); hold on

eval('fdm_3d_vargrid');
nErr = relerr([E([1 2],:) [-0.5; -0.5/4]]);
loglog(N,nErr(1:end-1),'b-', 'Linewidth',2, 'color', [0 0 1]); hold on
eval('fdm_3d_dat');
nErr = relerr([E([1 2],:) [-0.5; -0.5/4]]);
loglog(N,nErr(1:end-1),'b--', 'Linewidth',2, 'color', [0.5 0.5 1]); hold on

 %eval('harm_fdm_1d')
 %nErr = relerr([E([1 2],1:end) [1; 3]]);
 %loglog(N,nErr(1:end-1),'-', 'Linewidth',2,'color',[0.5   0.5 0.5]); hold on

axis([3e0 1e3 1e-3 1e1])
legend('Location','SouthWest','2D variable grid', '2D uniform grid', '3D variable grid', '3D uniform grid')
%legend('Location','SouthWest','1D variable grid', '1D uniform grid', '2D variable grid', '2D uniform grid', '3D variable grid', '3D uniform grid')
%legend('Location','SouthWest','2D variable grid', '2D uniform grid', '3D variable grid', '3D uniform grid')
grid on; set(gca,'minorgrid','none')

%ylabel('$$\frac{|\lambda_n \mbox{-}  \lambda|}{\lambda}$$','interpreter','latex','Rotation',0)
ylabel('$$\Big | \frac{\lambda_n \mbox{-}  \lambda}{\lambda} \Big |$$','interpreter','latex','Rotation',0)
xlabel('number of grid points per axis')

print -depsc ../../img/h_convergence_vargrid_new.eps
title('FDM Convergence Improvement with Variable Grids')

