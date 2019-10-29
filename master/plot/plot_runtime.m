figure(1); clf;  subplot(1.5,1.5,1);
figure(2); clf;  subplot(1.5,1.5,1);
figure(3); clf;  subplot(1.5,1.5,1);
figure(4); clf;  subplot(1.5,1.5,1);
path(path,'../dat')
 %eval('fdm_1d_dat')
 %figure(1); loglog(N,TE,'r','Linewidth',1.5); hold on
 %figure(2); loglog(N.^1,TE,'r','Linewidth',1.5); hold on
 %figure(3); loglog(N.^1,N.^1,'r','Linewidth',1.5); hold on
 %figure(4); loglog(N.^1,TE,'r','Linewidth',1.5); hold on
%eval('fdm_2d_dat')
eval('fdm_2d_vargrid_dat_new')
figure(1); loglog(N,TE,'g','Linewidth',1.5); hold on
figure(2); loglog(N.^2,TE,'g','Linewidth',1.5); hold on
figure(3); loglog(N.^1,N.^2,'g','Linewidth',1.5); hold on
figure(4); loglog(N.^2,TE,'g','Linewidth',1.5); hold on
eval('fdm_3d_dat')
figure(1); loglog(N,TE,'b','Linewidth',1.5); hold on
figure(2); loglog(N.^3,TE,'b','Linewidth',1.5); hold on
figure(3); loglog(N.^1,N.^3,'b','Linewidth',1.5); hold on
figure(4); loglog(N.^3,TE,'b','Linewidth',1.5); hold on
 %eval('fem_1d_dat')
 %figure(1); loglog(N,TE,'r--','Linewidth',1.5); hold on
 %figure(2); loglog(N.^1,TE,'r--','Linewidth',1.5); hold on
eval('fem_2d_dat')
%figure(1); loglog(N,TE,'g--','Linewidth',1.5); hold on
figure(2); loglog(N.^2,TE,'g--','Linewidth',1.5); hold on
eval('fem_3d_dat')
%figure(1); loglog(N,TE,'b--','Linewidth',1.5); hold on
figure(2); loglog(N.^3,TE,'b--','Linewidth',1.5); hold on

 %eval('fdm_1d_vargrid')
 %figure(4); loglog(N.^1,TE,'r--','Linewidth',1.5); hold on
eval('fdm_2d_vargrid')
figure(4); loglog(N.^2,TE,'g--','Linewidth',1.5); hold on
eval('fdm_3d_vargrid')
figure(4); loglog(N.^3,TE,'b--','Linewidth',1.5); hold on

figure(1)
 grid on;
 set(gca,'minorgrid','none');
 set(gca,'xtick',10.^(0:6));
 xlim([1 1e3]);
% xlim([1 1e6]);
 legend('location','southeast','2D FDM', '3D FDM', '2D FEM', '3D FEM');
% legend('location','southeast','1D FDM', '2D FDM', '3D FDM', '1D FEM', '2D FEM', '3D FEM');
 xlabel('number of gridpoints per axis');
 ylabel('time to compute first ten eigenvalues [s]');
 title('Eigensolver Runtime - Curse of Dimensionality')

figure(2)
 grid on;
 set(gca,'minorgrid','none');
 set(gca,'xtick',10.^(0:6));
 xlim([1 1e6]);
 legend('location','southeast','1D FDM', '2D FDM', '3D FDM', '1D FEM', '2D FEM', '3D FEM');
 xlabel('total number of gridpoints');
 ylabel('time to compute first ten eigenvalues [s]');
 title('Eigensolver Runtime - Curse of Dimensionality')

figure(3)
 grid on;
 set(gca,'minorgrid','none');
 set(gca,'xtick',10.^(0:6));
 xlim([1 1e6]);
 legend('location','southeast','1D', '2D', '3D');
 xlabel('gridpoints per axis n');
 ylabel('total number of grid points');
 title('Total Number of Grid Points - Curse of Dimensionality')

figure(1);
print -depsc ../../img/h_runtime_new.eps

