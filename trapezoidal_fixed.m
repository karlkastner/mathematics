% 2022-05-02 13:03:22.456110293 +0200
% interate in time using the trapezoildal scheme
function [t,yy,kk] = trapezoidal_fixed(fun,T,y0,opt)
	t    = T;
	nt   = length(t);
	ny = length(y0);
	jfun = @opt.Jacobian;
	tolgn = 0.001*sqrt(eps);
	kmax = ny;
	
	% initialize
	yy = zeros(nt,ny);
	yy(1,:) = y0;
	y = y0;
	df_dt = fun(t(1),y);
	I = speye(ny);
	kk=zeros(nt,1);
	for idx=2:nt
		dt = t(idx)-t(idx-1);

		% y_ = y + 1/2 dt (f(y_) + f(y))
		% initial guess
		y_ = y;
		% gauss newton iteration
		k = 0;
		while (true)
			k = k+1;
			% gradient
			df_dt_   = fun(t(idx),y_);
			J        = jfun(t(idx),y_);
			% residual
			r  = y_ - y - 0.5*dt*(df_dt_ + df_dt);
			dR_dy = (I - 0.5*dt*J);
			if (0)
				dy = dR_dy \ r;
			else
				dy = bicgstab(dR_dy,r,[],100);
			end
			y_ = y_ - dy;
			if (rms(dy) < tolgn)
				break;
			end 
			if (k>kmax)
				warning('no convergence');
				return;
			end
		end
		y = y_;
		df_dt = df_dt_;
		yy(idx,:) = y;
		kk(idx) = k;
	end % for idx
end

