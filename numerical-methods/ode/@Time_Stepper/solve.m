% Fri  7 Aug 17:46:01 +08 2020
% save set of odes by the adams-bashford multistep method
%
% solve ode and pdes by integrating over time
%
% function [to,fo] = solve(odefun,f0);
%
function [to,fo] = solve(obj,odefun,f0)

	nx = size(f0,1);

	% determine number of historic values held in temporary storage
	% TODO, the adams bashforth shemes need to be generalized for higher variable time steps
	switch (obj.scheme)
	case {'adams-bashforth'}
		nf  = 1;
		nd  = obj.norder;
		ddir = 0;
	case {'implicit'}
		nf = 1;
		nd = obj.norder;
		ddir = 0;
	case {'leapfrog'}
		nf = 2;
		nd = 1;
		ddir = 0;
	case {'leapfrog-trapezoidal'}
		nf = 3;
		nd = 1;
		ddir = 0;
	case {'upwind'}
		nf = 1;
		nd = 1;
		ddir = 1;
	otherwise
		error('here');
%		mode   = 'normal';
	end

	% allocate function value and derivative vector
	f     = zeros(nx,nf);
	df_dt = zeros(nx,nd);

	% slot indices for historic values held for the multi-step methods
	idd = 1:nd;
	idf = 1:nf;

	% initialize function
	f(:,idf(1))  = f0;

	% initialize derivative
	[df_dt(:,idd(1)),dt_max] = odefun(obj.Ti(1),f0,ddir);
	dt                       = obj.cfl*dt_max;
	dt_old = dt;

	% extrapolate history
	for idx=2:nf
		f(:,idf(idx)) = f(:,idf(1)) - (idx-1)*dt*df_dt(:,idd(1));
	end % for idx

	% initialize output function values
	fo      = [f0,zeros(nx,obj.n_realloc)];
	to      = [obj.Ti(1);zeros(obj.n_realloc,1)];
	odx = 1;
	t   = obj.Ti(1);
	% step in loop until final time is reached
	ti = tic();
	tstat = toc(ti)+10;
	while (t-obj.Ti(1)<(obj.Ti(2)-obj.Ti(1))*(1-10*eps))
		dt = obj.cfl*dt_max;

		% limit last time step
		dt = min(obj.Ti(2)-t,dt);
		if (toc(ti) > tstat)
			disp((t-obj.Ti(1))/(obj.Ti(end)-obj.Ti(1)));
			tstat = tstat+10;
		end	

		% step in time
		t = t+dt;
		switch (obj.scheme)
		case {'implicit'}
			% TODO df_dt does not need to be evaluated again later for the inmplicite solver 
			df_dt_p = df_dt;
			while (1)
				switch (obj.order)
				case {1} % euler
					fp = f + dt*df_dt_p;
				case {2} % trapezoidal
					fp = f + 0.5*dt*(df_dt_p + df_dt);
				otherwise
					error('not yet implemented');
				end % witch
				% derivative at next time-step
				df_dt_new = odefun(t+dt,fp,ddir);
				res       = df_dt_new - df_dt_p;
				if (rms(res) < obj.reltol*rms(df_dt_p))
					break;
				end
				df_dt_p   = df_dt_new;
			end % while (1)
			f = fp;
		case {'leapfrog-trapezoidal'}
			% the lft replace f2 by 0.5*(f1 + f3), or the equivalent
			% interpolate by for a variable time step
			fc      = f(:,idf(1)) ...
                                  + dt/(dt+dt_old)*(f(:,idf(3)) - f(:,idf(1)));
			% index 3 becomes index 1 after shifting
			f(:,idf(3)) = fc + 2*dt*df_dt;
			dt_old      = dt;
		case {'leapfrog'}
			% leap-frog
			% TODO interpolate
			f(:,idf())  = f(:,idf(2)) + 2*dt*df_dt(:,idd(1));

			% TODO robert-asselin-filter
		case {'mccormack'}
			% TODO implement
			% f = mccormack_step(tt,x,f,odefun,dt);
			error('not yet implemented');
		case {'adams-bashforth'}
			% TODO extend to variable time step
			switch (min(obj.order,tdx-1))
			case {1}
				% euler
				f = f + dt*df_dt(:,idd(1));
			case {2} 
				f = f + dt/2*(    3*df_dt(:,idd(1)) ...
	                                       -    df_dt(:,idd(2)) ...
	                                     );
			case {3}
				f = f + dt/12*(  23*df_dt(:,idd(1)) ...
	                                       - 16*df_dt(:,idd(2)) ...
					       +  5*df_dt(:,idd(3)) ...
				              );
			case {4}
				f =  f + dt/24*(   55*df_dt(:,idd(1)) ...
	                                         - 59*df_dt(:,idd(2)) ...
	                                         + 37*df_dt(:,idd(3)) ...
				                 -  9*df_dt(:,idd(4)) ...
						);
			otherwise
				error('not yet implmented')
			end % switch order
		case {'upwind'}
			f = f + dt*df_dt;
		otherwise
			error('not yet implemented');
		end % switch scheme

		% write state to output
		if (t > to(odx)+obj.dto)
			odx=odx+1;
			% reallocate
			if (size(fo,2)<odx)
				fo(:,end+(1:obj.n_realloc)) = 0;
				to(end+(1:obj.n_realloc))   = 0;
			end
			fo(:,odx) = f(:,idf(1));
			to(odx)   = t;
		end

		% shift indices, instead of the derivative vectors
		idd = circshift(idd,1);
		idf = circshift(idf,1);
	
		% evaluate derivative at current state
		[df_dt(:,idd(1)), dt_max] = odefun(t,f(:,idf(1)),ddir);

	end % for tdx

	% write final state
	odx       = odx+1;
	fo(:,odx) = f(:,idf(1));
	to(:,odx) = t;

	% truncate
	fo = fo(:,1:odx);
	to = to(1:odx);
end % solve

