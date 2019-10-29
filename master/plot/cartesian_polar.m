 l=linspace(0,2*pi,20);
 hold on;
 plot3(sin(l), cos(l), zeros(size(l)), 'k', 'linewidth', 2), quiver3(0, 0, 0, [sqrt(1/3)],[sqrt(1/3)], [sqrt(1/3)]);
 quiver3(0, 0, 0, [sqrt(1/2)],[sqrt(1/2)], 0);
 quiver3(0, 0, 0, 1, 0, 0);
 l2 =  linspace(0,2*pi/8,20);
 plot3(0.5*sin(l2+pi/4), 0.5*cos(l2+pi/4), zeros(size(l2)), 'k', 'linewidth', 2), quiver3(0, 0, 0, [sqrt(1/3)],[sqrt(1/3)], [sqrt(1/3)]);
 plot3(0.5*sin(l2), 0.5*sin(l2), 0.5*(cos(l2-pi/2)), 'k', 'linewidth', 2)
 quiver3(0, 0, 0, [sqrt(1/3)],[sqrt(1/3)], [sqrt(1/3)]);

 set(gca,'visible','off')

 [X Y Z] = sphere(20);
 h=surf(X,Y,Z);
 set(h,'facecolor',[0.5 0.5 0.5],'facealpha',0.5,'edgecolor','none');
 set(h,'FaceLighting','phong');
 axis equal;
 light('Position',[1 0 0],'Style','infinite');
                                       

 ms = 16;
 n = 20;
 fs = 16;
 s = 1;
 l = linspace(0,s,n);
 o = s*ones(n,1)
 subplot(2,2,1)
 plot(o,l,'k--'); hold on
 plot(l,o,'k--')
 plot(s,s,'k.','markersize',ms)
 axis(s*[0 1.5, 0 1.5])
 text(1.1*s,0.50*s,'y','FontSize',fs)
 text(0.50*s,1.1*s,'x','FontSize',fs)
 set(gca,'xtick',[]);
 set(gca,'ytick',[]);

 subplot(2,2,2);
 s2 = 0.5*s;
 plot(l,l,'k--'); hold on
 l2 = linspace(0,pi/4,n);
 plot(s2*cos(l2),s2*sin(l2),'k--')
 plot(s,s,'k.','Markersize',ms);
 axis(s*[0 1.5, 0 1.5])
 text(0.9*s/2,1.2*s/2,'r','FontSize',fs)
 text(2/8,2/16,'$\phi$','interpreter','latex','FontSize',fs)
 set(gca,'xtick',[]);
 set(gca,'ytick',[]);
 print -deps cartesian-polar.eps
 system('epstopdf cartesian-polar.eps')
