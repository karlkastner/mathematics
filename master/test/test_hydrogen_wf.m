function test_hydrogen_wf

opengl neverselect

clf
L = 10;
nr = 50;
nt = nr;
np = nr;%/2;
r     = linspace(-L,L,nr); 
theta = linspace(0,pi,nt);
phi   = linspace(0,2*pi,np);

n_max = 4;
p_= 1;
for n=1:n_max
 for l=0:n-1
  for m=0:l
	[n l m]
	[Psi RR TT PP X Y Z] = hydrogen_wf(n,l,m,r,theta,phi);
	
	subplot(n_max,n_max+1,p_);
	p_ = p_+1;
	lim = max(max(max(Psi)));
	%p = patch(isosurface(XX,YY,ZZ,Psi,lim));
	p0 = patch(isosurface(X,Y,Z,Psi, 0.5*lim));
	set(p0,'FaceColor','yellow','EdgeColor','none');
	p1 = patch(isosurface(X,Y,Z,Psi, 0.05*lim));
	set(p1,'FaceColor','green','EdgeColor','none');
	p2 = patch(isosurface(X,Y,Z,Psi, 0.005*lim));
	set(p2,'FaceColor','blue','EdgeColor','none');
	p3 = patch(isosurface(X,Y,Z,Psi, 0.0005*lim));
	set(p3,'FaceColor','red','EdgeColor','none');

	daspect([1,1,1])
	view(3);
	camlight
	lighting gouraud
%	grid on
	view([90 90])
	title(num2str([n l m]))
	axis([-25 25 -25 25 -25 25])
	set(gca,'XTickLabel',[])
	set(gca,'YTickLabel',[])
	set(gca,'ZTickLabel',[])

	pause(1)
	
  end % for m
 end % for n
end % for l

print -depsc ../img/hydrogen_wf_3d.eps

end

