% Wed 29 Mar 13:43:02 CEST 2017
% Karl Kastner, Berlin
%
%% approximate range of a continous Fourier series with 2 components
%% range(y) = max(y) - min(y)
% range
% midrange
% phi0
function [range, midrange, y0, phi0] = fourier_range(a,phi)
	a0 = a(1,:);
	a  = a(2:3,:);

	phi_     = phi;
	dp       = phi(2,:)-2*phi(1,:);
	if (~issym(dp))
	dp       = wrapTo2Pi(dp);
	end
	phi(1,:) = 0;
	phi(2,:) = dp;

	phi0 = [2*pi-(a(1,:).*phi(1,:) + 2*a(2,:).*phi(2,:))./(a(1,:)+4*a(2,:))
	        2*pi-(a(1,:).*(phi(1,:)+pi) + a(2,:).*(2*phi(2,:)+2*pi))./(a(1,:)+4*a(2,:))];

	if (~issym(a))
	phi0_ = [-(a(1,:).*phi(1,:) + 2*a(2,:).*(-2*pi+phi(2,:)))./(a(1,:)+4*a(2,:))
	          2*pi-(a(1,:).*(phi(1,:)+pi) + a(2,:).*(2*phi(2,:)+2*pi))./(a(1,:)+4*a(2,:))];

	fdx = wrapTo2Pi(phi(2,:) - 2*phi(1,:)) > pi;
	phi0(1,fdx) = phi0_(1,fdx);
	end

	if (~issym(a))
	phi0 = bsxfun(@plus,phi0,-phi_(1,:));
	else
		phi0 = phi0-phi_;
	end

	if (~issym(a))
	phi0 = wrapTo2Pi(phi0);
	end
	
	phi=phi_;
	y0 = [a(1,:).*cos(phi0(1,:)+phi(1,:))+a(2,:).*cos(2*phi0(1,:)+phi(2,:));
	      a(1,:).*cos(phi0(2,:)+phi(1,:))+a(2,:).*cos(2*phi0(2,:)+phi(2,:)) ];

	range    = -(y0(1,:)-y0(2,:));
	midrange = a0 + 1/2*sum(y0); 
end

