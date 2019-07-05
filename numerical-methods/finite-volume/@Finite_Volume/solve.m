% Thu Apr 28 04:43:13 MSD 2011
% Sun Jun 30 18:00:09 UTC 2013
% 2013/04/14 23:45
% Mi 3. Feb 12:27:33 CET 2016
% Karl Kästner, KTH Stockholm
%
%% solve the the PDE by successively stepping in time
%% this is a trivial implmentation with constant step length
%% severity of diffusive error depends on dt/dx-ratio
%% stability depends on wave height
%
% adaptive time step pde solver
function [T, Q] = solver(obj, Ti) %, b0fun)
	n = length(obj.x);

	startup = obj.startup;
%	abstol  = 1e-3;
%	abstol    = 1e-7;
%	reltol    = 1e-5;
	idstartup = 0;

	% initial condition
	q      = obj.icfun(obj.x);

	t       = Ti(1);
	tout    = Ti(1);
	odx     = 1;
	tdx     = 1;
	dt_new  = obj.cflscale*obj.pde.dt_cfl(q,obj.dx);
	timer   = tic();
	t_progress = toc(timer);


	% start values
	Q(:,1) = q;
	T      = t;
	while (t<Ti(end))
		if (any(isnan(q)) || any(q(1:end/2)<0)) %dt_new))
			error('NaN or negative water depth occured, at step %d cannot continue\n',tdx);
			break;
		end

		tdx=tdx+1;
		% limit step to end of time series
		dt_new = min(dt_new,(Ti(end)-t)*(1+sqrt(eps)));

		% try to step
		while (1)
			% shrink towards zero
			dt = obj.cflscale*dt_new;

			dt = min(dt,tout-t+obj.dt_out);

		        q_new = obj.step(t, q, dt);


			% check cfl
			dt_new = obj.pde.dt_cfl(q_new,obj.dx);
 			if (any(isnan(q_new)) || any(q_new(1:end/2)<0))
				error('Solution is nan');
				q = q_new;
				break;
			end
			% check
			if (dt < dt_new & startup)
				idstartup = idstartup+1;
				%delta = max(abs(q-q_new))/dt;
				delta_new = q_new-q;
				if (idstartup > 1)
					% q_new = linacc(q_new,q,delta_new,delta);  
				end

				%delt = max(delta);
				% p         = delta/delta_old;
				% delta_inf = delta*delta_old/(delta_old-delta);
				% delta_old = delta;
				%if (delta_inf >= 0 && delta_inf < obj.icabstol)
				% all(delta < abstol | delta < reltol*abs(q)))
				delta = max(abs(delta_new));
				if (all(delta < obj.icabstol | delta < obj.icreltol*abs(q)))
					startup = false;
					Q(:,1) = q;
					printf('Startup required %d iterations\n',idstartup);
%			figure(1)
%			clf
%			subplot(2,1,1)
%			plot([q(1:end/2) q_new(1:end/2)])	
%			subplot(2,1,2)
%			plot([q(end/2+1:end) q_new(end/2+1:end)])	
%				pause

				else
					delta = q_new - q; 
					q=q_new;
					continue;
				end
			end
			if (dt < dt_new)
				% accept and step in time
				q  = q_new;
				t  = t + dt;
				break;
%			end
%			if (isnan(dt))
%				q  = q_new;
%				t  = t + dt;
%				break;
			end
		end % while 1

		% store result
		if (t >= tout+obj.dt_out)
			odx      = odx+1;
			tout     = t;
			% reallocate memory
			if (obj.dt_out>0)
				nt       = odx+ceil((Ti(end)-t)/obj.dt_out);
			else
				nt       = odx+ceil((Ti(end)-t)/dt);
			end
			if (nt > size(Q,2))
				Q(1,nt) = NaN;
				T(nt,1) = NaN;
			end
			T(odx)   = t;
			Q(:,odx) = q;
		end

		% output
		% obj.measure_points(T(idx),A,Q);
		% obj.measure_map(T(idx),A,Q);
		t_real = toc(timer);	
		if (t_real > t_progress + obj.dt_progress)
			printf('Progress %2.1f%% %2.1fs\n',100*(t-Ti(1))/(Ti(2)-Ti(1)),t_real);
			t_progress = t_real;
		end
        end % while t < T(end)

	% truncate, if more memory was allocated then required
	Q = Q(:,1:odx);
	T = T(1:odx);

	% depth to surface elevation
	%Q(1:n,:) = Q(1:n,:) + repmat(zb,1,length(T));
end % solver

% linear acceleration
function x0 = linacc(x1,x2,y1,y2)
	% y = y1+(x-x1)*(y2-y1)/(x2-x1)
%	y1 = y1'*y1;
%	y2 = y2'*y2;
	fdx = abs(y2-y1) > sqrt(eps);
%	x0 = 
	%dydx = (y2-y1)./(x2-x1);
	
	p = 0.01;
	
	dxdy = (x2-x1)./(y2-y1);
	x0 = p*x1 + (1-p)*(x1 - y1.*dxdy);
	x0(~fdx) = x1(~fdx);
	%x0 = x1 - y1./dydx;
%	x0 = x1;
%	x0(fdx) = x1(fdx).*(y2(fdx)-y1(fdx))-y1(fdx).*(x2(fdx)-x1(fdx))./(y2(fdx)-y1(fdx));
%	x0 = x1.*(y2-y1)-y1.*(x2-x1)./(y2-y1);
end

