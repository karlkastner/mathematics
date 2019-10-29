% Mon Feb 20 19:03:54 MSK 2012
% Karl Kästner, Berlin

 clf
 subplot(2,2,1)
 plot([-0.5 0 0.5 1 1.5 ], [1 0 1 0 1], 'k--','Linewidth',2)
 hold on
 plot([0 0.5 1 ], [0 1 0], 'k-','Linewidth',2)
 plot([-0.5 0 0.5 1 1.5], [0 1 0 1 0], 'k--','Linewidth',2)
 text(0.1,0.9,'\fontsize{16}\phi_{i-1}')
 text(0.6,0.9,'\fontsize{16}\phi_{i}')
 text(1.1,0.9,'\fontsize{16}\phi_{i+1}')

 text(0.075,0.05,'\fontsize{16}x_{i-1}')
 text(0.575,0.05,'\fontsize{16}x_{i}')
 text(1.075,0.05,'\fontsize{16}x_{i+1}')
 
 text(0.0, 1.3,'\fontsize{16}u_{i-1}')
 text(0.5, 1.6,'\fontsize{16}u_{i}')
 text(1.0, 1.70,'\fontsize{16}u_{i+1}')

 plot([-0.5 0 0.5 1 1.5], [1 1.125 1.375 1.875 1.625], 'k.-','Linewidth',2,'Markersize',16)

 axis([-0.1 1.3 0 2]); 
 set(gca,'xticklabel',[]);
%xlabel('Honk');ylabel('Molch');title('Kalle');
 print -deps ../img/hat_function.eps

