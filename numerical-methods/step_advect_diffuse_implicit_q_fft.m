% Thu 29 Feb 15:14:46 CET 2024
% TODO
function zout = step_advect_diffuse_implicit_q_fft(dt,dx,a,e,z0,q,isreal_)
	if (isvector(z0))
	n = length(z0);
	a = a./dx;
	e = e./dx.^2;
	
  % frequencies
 	% o = 2*pi*[(0:n/2),(1:n/2-1)]';                                                     
	o = 2*pi*[(0:n/2-1),-((n/2):-1:1)]';
%          o =2*pi*[0:n/2, -(n/2-1:-1:1)]
	% eigenvalues
	ee = (   e(1)*2*(cos(o/n)-1) ...
	       + a(1)*(+1 - cos(o/n) + 1i*sin(o/n)).*(a(1)<0) ...
	       + a(1)*(-1 + cos(o/n) + 1i*sin(o/n)).*(a(1)>0) ...
	     );
 	           % (1 - cos(d/n) - 1i*sin(-d/n) )
%  (I - a*D2)z = z0
%  V*V'(I - q*a*D2)*V*V'*z = V*V'(I + q*a*D2)*V*V'*z0
%    V'(I - a*D2)*V*V'*z =   (I + q*a*E)*V'*z0
%    V'(I - a*D2)*V*tz   =   (I+q*a*E)*tz0
%      (I - a*V'*D2*V)*tz   =  (I+q*a*e)*tz0
%      (I - a*E)^-1*tz   =   tz0
%      tz = tz0/(1-a*E)
%      z  = V*(tz0/(1-a*E))
	tz0 = (fft(z0)/sqrt(n));
   	% evolve in time
   	tz = tz0.*(1+(1-q)*dt*ee)./(1-q*dt*ee);
   	% retransform
   	zout = ifft(sqrt(n)*tz);

	else
		a = cvec(a)./cvec(dx);
		e = cvec(e)./cvec(dx).^2;

		n = size(z0);
		%o1 = 2*pi*[(0:n(1)/2-1),-((n(1)/2):-1:1)]'/n(1);
		%o2 = 2*pi*[(0:n(2)/2-1),-((n(2)/2):-1:1)]/n(2);
		[fx,fy] = fourier_axis_2d(n,n);
		o1 = 2*pi*fx;
		o2 = 2*pi*fy';

		p = 1;
		ee = (  2*e(1)*(cos(o1)-1) ...
	              + 2*e(2)*(cos(o2)-1) ...
		      + a(1)*( 1 - cos(o1) + 1i*sin(o1)).*(a(1)<0) ...
		      + a(1)*(-1 + cos(o1) + 1i*sin(o1)).*(a(1)>0) ...
		      + a(2)*( 1 - cos(o2) + 1i*sin(o2)).*(a(2)<0) ...
		      + a(2)*(-1 + cos(o2) + 1i*sin(o2)).*(a(2)>0) ...
		     );

		% note, the normalization cancels
		% tz0 = (fft2(z0)/sqrt(n(1)*n(2)));
		tz0 = (fft2(z0)); % 1/sqrt(n(1)*n(2)));
   		tz = tz0.*(1+((1-q)*dt)*ee)./(1-(q*dt)*ee);
		zout  = ifft2(tz); % *sqrt(n(1)*n(2)));
	end                            
	if (isreal_)
		zout = real(zout);
	end
end                                                       

