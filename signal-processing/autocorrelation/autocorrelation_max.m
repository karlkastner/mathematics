% Sun 15 Jan 21:48:27 CET 2023
% radial density along main axis and ring passing through first side lobe
% for isotropic patterns

% TODO only when dx=dy, Ly can differ from Lx
function [Rr,r,Rt,t,Rc,xc,yc,out] = autocorrelation_max(R,L)
	
	n = size(R);
	
	% move origin to centre
	R = fftshift(R);
	R0 = R;
	
	% strip symmetric upper half
	% if odd, then on (floor(n/2)+1) will be the centre
	R(1:floor(n(1)/2),:) = [];
	
	% mask the main lobe
	% to right
	idx = ceil(n(2)/2);
	while (true)
		jdx = 1;
		if (R(jdx,idx) <= 0)
			break;
		end
		while (true)
			if (R(jdx,idx) <= 0)
				break;
			end
			R(jdx,idx) = 0;
			jdx = jdx+1;
		end
		idx = idx+1;
	end
	
	% to left
	idx = ceil(n(2)/2)-1;
	while (true)
		jdx = 1;
		if (R(jdx,idx) <= 0)
			break;
		end
		while (true)
			[jdx,idx]
			if (R(jdx,idx) <= 0)
				break;
			end
			R(jdx,idx) = 0;
			jdx = jdx+1;
		end
		idx = idx-1;
	end
	
	% find the maximum
	[Rc,ijc] = max(R,[],'all');
	
	% transform index to subscripts
	[ic,jc] = ind2sub(size(R),ijc);
	jc = jc - ceil(n(2)/2);
	xc = jc*L(2)/n(2);
	yc = ic*L(1)/n(1);
	h = hypot(ic,jc);
	
	% right quarter
	% if (jc>n(2)/2)
	%	jc = n(2)-jc-1;
	% end
	
	% transform to central coordinates
	%R = fftshift(R);
	%ic = ic+n(1)/2;
	%jc = jc+n(2)/2;
	
	% extract autocorrelation along max
	c = ic/h;
	s = jc/h;
	A = [c,s;
	       -s,c];
	t0 = atan2(c,s);
	nmax = round(min(n)/2);
	ri = (0:nmax-1);
	r  = ri*L(1)/n(1);
	ij = A*flipud([ri;
	        zeros(1,nmax)]) + [floor(n(2)/2);floor(n(1)/2)]+1;
	%ij(1,:)+ceil(n(1)/2),ij(2,:)+ceil(n(2)/2));
	Rr = improfile(R0,ij(1,:),ij(2,:),'bilinear');
	
	nt  = round(2*pi*h);
	t   = 2*pi*(0:nt-1)/nt;
	i   = floor(n(1)/2)+1 + h*sin(t+t0); 
	j   = floor(n(2)/2)+1 + h*cos(t+t0);
	
	Rt  = improfile(R0,j',i','bilinear');
	
	out.ic = ic;
	out.jc = jc;
	out.rij = ij;
	out.tij = [i;j];

end

