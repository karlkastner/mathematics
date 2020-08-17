% Fri  7 Aug 17:46:01 +08 2020
% save set of odes by the adams-bashford multistep method
%
% function [tt,ff] = ode_adams_bashforth(odefun,f0,T,nt,norder,ks);
% TODO, check stability
% note : this method cannot vary the time step
function [tt,ff] = ode_adams_bashforth(odefun,f0,T,nt,norder,ks,steptol);
	correction = false;
	reltol = 1e-3;
	lambda = 0.25;

	nx = size(f0,1);

	switch (norder)
	case {'leap-frog','leapfrog','leapfrog-trapezoidal','lft'}
		mode = 'leapfrog';
		norder = 1;
	case {'upwind'}
		norder = 1;
		mode = 'upwind';
	otherwise
		mode = 'normal';
	end

	% initialize derivative
	df_dt = zeros(nx,norder);

	% slots for last nprder-derivatives held for the multi-step
	id = 1:norder;


	f  = f0;
	% initialize
	if (~isempty(nt))
		dt = (T(end)-T(1))/nt;
		tt = T(1) + (T(end)-T(1)).*(0:nt)'./nt;
		nt_ = nt;
	end
	[df_dt(:,1),dtmax] = odefun(T(1),f0);
	if (isempty(nt))
		%dt = steptol/max(abs(df_dt));
		dt = steptol*dtmax;
	end
	f1 = f;
	f2 = f-dt*df_dt;
	f3 = f-2*dt*df_dt;
	%for tdx=2:nt+1
	% initialize output function values
	nt_ = 100;
	ff  = zeros(nx,ceil(nt_/ks));
	kdx = 1;
	ff(:,1) = f;
	tdx = 1;
	t = T(1);
	tt = zeros(nt_,1);
	pt = 0.5;
	while (t<T(2)*(1-10*eps))
		tdx = tdx+1;
		if (isempty(nt))
			dt = min(steptol*dtmax,2*dt);
			%dt = pt*steptol/max(abs(df_dt)) + (1-pt)*dt;
			dt = min(T(2)-t,dt);
		end
		disp(tdx)
%		dt
%		t/T(2)
		% step in time
		t = t+dt;
		switch (mode)
		case {'leapfrog-trapezoidal','lft'}
			f  = 0.5*(f1+f3) + 2*dt*df_dt(:,id(1));
			f3 = f2;
			f2 = f1;
			f1 = f;
		case {'leapfrog'}
			% leap-frog
			f  = f2 + 2*dt*df_dt(:,id(1));
			% robert-asselin-filter
			a  = 0.01;
			f1 = f1 + a*(f - 2*f1 + f2);
			f2 = f1;
			f1 = f;
		otherwise
		switch (min(norder,tdx-1))
		case {1}
			% euler
			if (0)
			f = mccormack_step(tt,x,f,odefun,dt);
			else
			fp = f + dt*df_dt(:,id(1));
			if (correction)
				df_dt_p_old = df_dt(:,id(1));
				while (1)
					df_dt_p_new = odefun(tt+dt,fp);
					res = df_dt_p_new - df_dt_p_old;
					df_dt_p = lambda*df_dt_p_new + (1-lambda)*df_dt_p_old;
%)./norm(df_dt_p)			df_dt_p_ = df_dt_p;
					% semi-implicit step
					fp = f + dt*df_dt_p;
					e = rms(res)./rms(df_dt_p_new);
					e*lambda
					if (e<lambda*reltol)
						break;
					end
					df_dt_p_old = df_dt_p_new;
				end
			end
				f = fp;
			end
			
		case {2} 
			f = f + dt/2*(3*df_dt(:,id(1)) - df_dt(:,id(2)));
		case {3}
			f = f + dt/12*(23*df_dt(:,id(1)) - 16*df_dt(:,id(2)) + 5*df_dt(:,id(3)));
		case {4}
			f = f + dt/24*(55*df_dt(:,id(1)) - 59*df_dt(:,id(2)) + 37*df_dt(:,id(3)) - 9*df_dt(:,id(4)));
		otherwise
			error('not yet implmented')
		end
		end

		% save
		if (~isempty(nt))
		if (tdx == (ks*kdx+1))
			kdx = kdx+1;
			ff(:,kdx) = f;
		end
		else
			kdx=kdx+1;
			if (size(ff,2)<kdx)
				ff(:,end+(1:nt_)) = 0;
				tt(end+(1:nt_)) = 0;
			end
			ff(:,kdx) = f;
			tt(kdx)   = t;
		end

		% shift indices, instead of the derivative vectors
		id = circshift(id,1);
	
		% evaluate derivative
		[df_dt_, dtmax] = odefun(t,f);
		p_ = 1;
		df_dt(:,id(1)) = p_*df_dt_ + (1-p_)*df_dt(:,id(1)); 
%		df_dt(id)
%		2*pi*cos(2*pi*tt(tdx:-1:max(tdx-3,1)))
%pause
	end % for tdx
	kdx = kdx+1;
	ff(:,kdx) = f;
	tt(:,kdx) = t;
	% truncate
	ff = ff(:,1:kdx);
	tt = tt(1:kdx);
end

