clf
m = 6;
for idx=2:m
    subplot(1,1.5*(m-1),1.5*(idx-1))
    n = 1;
    for jdx=1:idx
	plot(((1:n)-idx)/(idx-1)+ 0.5*(idx-n)/(idx-1),(idx-jdx)/(idx-1),'k.','MarkerSize',10); hold on
	n = n+1;
    end
%   text(0.5*n/idx,-1/idx,'Honk');
%     set(gca,'xtick',[]);
%     set(gca,'ytick',[]);
%     axis([-0.1 1.1 -0.1 1.1]);
     axis equal
     axis tight;
%     axis tight
%    set(gca,'xcolor','w')
%    set(gca,'ycolor','w')
%    xlabel(n*(n-1)/2);
    set(gca,'visible','off')
    text(-0.6,-0.2,num2str(n*(n-1)/2));

%    axis equal;
 % box off;
%    axis off
%    box on
%    axis([-0.6 0.6 -0.1 1.1])
    %set(gca,'Linewidth',0)
%    set(gca,'visible','off')
    
%   title(n*(n-1)/2);
%    text(0.5,-0.1, num2str((n*(n-1)/2)));
end
 
print -deps triangle-dof.eps
