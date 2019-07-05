% Wed 16 May 12:05:52 CEST 2018
% Karl Kastner, Berlin
%
%% fourier series coefficient of a ramp
%
function ab = fourier_coefficient_ramp3(n,L,x0,x1)
	afun =     @(L,n,x0,x1)(-L.*sin((pi.*n.*x0.*2.0)./L)+L.*sin((pi.*n.*x1.*2.0)./L)-L.*sin((pi.*n.*(L.*2.0-x0).*2.0)./L)+L.*sin((pi.*n.*(L.*2.0-x1).*2.0)./L)+L.*pi.*n.*cos((pi.*n.*x0.*2.0)./L).*4.0-L.*pi.*n.*cos((pi.*n.*x1.*2.0)./L).*4.0+pi.^2.*n.^2.*x0.*sin(pi.*n.*2.0).*8.0-pi.^2.*n.^2.*x1.*sin(pi.*n.*2.0).*8.0)./(pi.*n.*(x0-x1).*(cos(pi.*n.*4.0)+pi.^2.*n.^2.*8.0-1.0));

    bfun = @(L,n,x0,x1)(-L.*cos((pi.*n.*x0.*2.0)./L)+L.*cos((pi.*n.*x1.*2.0)./L)+L.*cos((pi.*n.*(L.*2.0-x0).*2.0)./L)-L.*cos((pi.*n.*(L.*2.0-x1).*2.0)./L)+L.*pi.*n.*sin((pi.*n.*x0.*2.0)./L).*4.0-L.*pi.*n.*sin((pi.*n.*x1.*2.0)./L).*4.0-pi.*n.*x0.*sin(pi.*n.*2.0).*4.0+pi.*n.*x1.*sin(pi.*n.*2.0).*4.0-pi.^2.*n.^2.*x0.*cos(pi.*n.*2.0).*8.0+pi.^2.*n.^2.*x1.*cos(pi.*n.*2.0).*8.0)./(pi.*n.*(x0-x1).*(cos(pi.*n.*4.0)+pi.^2.*n.^2.*8.0-1.0));

	ab(:,1) = afun(L,n,x0,x1);
	ab(:,2) = bfun(L,n,x0,x1);

	fdx = (0 == n);
	ab(fdx,1) = (1 -x0/L - 1/2*(x1-x0)/L - 0*(L-x1)/L);
	ab(fdx,2) = 0;

