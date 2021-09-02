% Wed 28 Jul 19:37:55 CEST 2021
function [r,Rmax,xmax] = acf_decay_rate(x,R)
	R = R(1:end/2);
	x = x(1:end/2);
	fdx = find(R<0,1,'first');
	% first lobe
	fdx1 = find(R(fdx+1:end)>0,1,'first')+fdx;
	fdx2 = find(R(fdx1+1:end)<0,1,'first')+fdx1;
	[Rmax,maxdx] = max(R(fdx1+1:fdx2-1));
	xmax = x(maxdx+fdx1);
	r = -log(Rmax)/xmax; 
end

