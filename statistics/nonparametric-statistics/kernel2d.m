% 2016-03-01 19:39:31.629833650 +0100 kernel2d.m
% Karl Kastner, Berlin
%
%% kernel density estimate in two dimensions
%function [Z, X, Y] = kernel2d(X,Y,xi,yi,sx,sy,rho,n1,n2,fun)
function [Z, X, Y] = kernel2d(X,Y,xi,yi,sx,sy,rho,n1,n2,fun)
	n    = length(xi);

	if (nargin() < 10)
		fun = @normpdf2;
	end

	% output coordinates
	if (isempty(X))
		if (nargin() < 8 || isempty(n1))
			n1   = round(sqrt(n));
		end
		xmin = min(xi);
		xmax = max(xi);
		dx   = xmax-xmin;
		xmin = xmin-0.5*dx;
		xmax = xmax+0.5*dx;
		X = linspace(xmin,xmax,n1)';
	else
		n1 = length(X);
	end
	if (isempty(Y))
		 if (nargin() < 9 || isempty(n2))
			n2   = n1;
		 end
		ymin = min(yi);
		ymax = max(yi);
		dy = ymax-ymin;
		ymin = ymin-0.5*dy;
		ymax = ymax+0.5*dy;
		Y = linspace(ymin,ymax,n2);
	else
		n2 = length(Y);
	end
	if (isvector(X))	
		X = cvec(X);
		X = repmat(X,1,n2);
	end
	if (isvector(Y))
		Y=rvec(Y);
		Y = repmat(Y,n1,1);
	end

	xyi = [cvec(xi) cvec(yi)];
	XY = [flat(X) flat(Y)];
	
	nxi = length(xi);
	% kernel density scale and covariance
	if (nargin() < 5 || isempty(sx))
		S2    = cov([xi,yi])/n;
		[v e] = eig(S2);
		S     = v*(8*sqrt(e)/sqrt(nxi))*v';
		S     = repmat(S,[1 1 nxi]);
%		S(1,1,1:nxi) = S(1,1);;
%		S(1,2,1:nxi) = rho*sx.*sy;
%		S(2,2,1:nxi) = sy.^p;
%		S(2,1,1:nxi) = rho*sx.*sy;
	else
		% todo power2  depends on kernel fun
		p = 1;
		S(1,1,1:nxi) = sx.^p;
		S(1,2,1:nxi) = rho*sx.*sy;
		S(2,2,1:nxi) = sy.^p;
		S(2,1,1:nxi) = rho*sx.*sy;
	end
	%Si = inv(S);

	% convolve
if (0)
	Z = zeros(n1*n2,1);
	for idx=1:length(xi)
		%Sidxdy = Si*[dx dy]';
		%z = dx* dy]'*Si*[dx dy];
		%z = 1./sqrt(2*pi*dS)*exp(-0.5*(dxy.*(Si*dxy')));
		z = fun(XY,xyi(idx,:),S);
		%z = normpdf([dx dy]*(Si*[dx dy]'));
		Z = Z + z; 
	end
	% normalise
	Z = 1/n*Z;
else
	scale = ones(length(xyi),1);
	for jdx=1:1
	Z = zeros(n1*n2,1);
	for idx=1:length(xi)
		z = fun(XY,xyi(idx,:),squeeze(S(:,:,idx))*scale(idx));
		%z = fun(XY,xyi(idx,:),S*scale(idx));
		Z = Z + z; 
	end
	% normalise
	Z = 1/n*Z;
	zi   = interp2(XY(1:n1,1),XY(1:n1:end,2),reshape(Z,n1,n2),xyi(:,1),xyi(:,2),'linear',min(Z(:)));
	p = 0.0;
	detS = 1; %det(S);
%	scale = normpdf(0,0,sqrt(detS))^2./(zi+sqrt(eps));
	scale = normpdf(0,0,S(1,1))*normpdf(0,0,S(2,2)) ./zi;
	%scale = (2*pi*sqrt(detS))*normpdf(0)^2./(zi+sqrt(eps));
	% normalise
	%min(scale)
	%scale = scale/min(scale);
	scale = scale/median(scale);
	scale = (1-p) + p*scale;
	end % for jdx
end


	X = reshape(X,n1,n2);
	Y = reshape(Y,n1,n2);
	Z = reshape(Z,n1,n2);
end

