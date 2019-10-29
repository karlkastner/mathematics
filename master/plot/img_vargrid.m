n = 12;
L = 12;
%X=((-n+1-1:n+1)-1/2)';
X=((-n+1:n)-1/2)';
e=1/pi;
Y = (X(end)/(e*(exp(e*X(end))-1))) * e * sign(X).*(exp(e*abs(X)) - 1);
%X = X(2:end-1);
%Y = Y(2:end-1);

clf
subplot(2,3,2)
line(L*[-ones(1,2*n); ones(1,2*n)],[Y'; Y'],'color',[0 0 0])
line([Y'; Y'], L*[-ones(1,2*n); ones(1,2*n)],'color',[0 0 0])
set(gca,'xtick',[])
set(gca,'ytick',[])
axis square
axis tight
axis([Y(1) Y(end) Y(1) Y(end)])
set(gca,'box','on')
subplot(2,3,1)
line(L*[-ones(1,2*n); ones(1,2*n)],[X'; X'],'color',[0 0 0])
line([X'; X'], L*[-ones(1,2*n); ones(1,2*n)],'color',[0 0 0])
set(gca,'box','on')
set(gca,'xtick',[])
set(gca,'ytick',[])
axis square
axis tight
axis([X(1) X(end) X(1) X(end)])

print -deps ../../img_vargrid.eps

