% Fri 14 Aug 16:01:18 +08 2020
function [A,b] = bvp1_row(cc,ll,cdx,dx,xi,bcfun)
	nxc = size(cc,1);
	% root for homogeneous part
	r  = ll(:,1,cdx-1);
	% discharge condition
	% set discharge for channels, where given
	[v1] = bcfun(xi(1),[],cdx);
	rid = 1;
	if (~isempty(v1))
		% set discharge according to boundary condition
		%A(2*nxc+1,2*nxc+1) = 1;
		A(1,2*nxc+1) = 1;
		% v1 == Q0
		b(2*nxc+1)         = -v1;
		%b(2*nxc+1)         = -v1;
		b(1)         = -v1;
	else
		[v2] = bcfun(xi(2),[],cdx);
		if (~isempty(v2))
			% set discharge according to boundary condition
				%A(2*nxc+1,2*nxc+1) = 1;
			A(1,2*nxc+1) = 1;
			% v1 == Q0
				%b(2*nxc+1)         = +v2;			
			b(1)         = +v2;			
		else
			% cd |Q0|Q0 + sum g hc^3 dz/dx = -1/2 cd |Q1|^2
			for jd=1:ncx
				% g*h*dz/dx = g*h*d/dx(a*exp(r*x)) = g*h*a*r
				% TODO, degenerated?
				%A(2*nxc+1,jd) = g*Ac(jd).*r(jd);
				A(rid,jd) = cc(jd,1).*r(jd);
				% Q0
				%A(2*nxc+1,2*nxc+1) = A(2*nxc+1,2*nxc+1) - cdc(jd)./(wc(jd).*hc(jd).^2).*abs(Q0);
				A(rid,2*nxc+1) = A(2*nxc+1,2*nxc+1) + cc(jd,2);
				b(rid)         = b(2*nxc+1) + cc(jd,3);
				%+ cdc(jd).*0.5.*abs(Q1(jd)))./(wc(jd).*hc(jd).^2);
			end % for
		end % if 
	end % if 
end % bvp1c_row


