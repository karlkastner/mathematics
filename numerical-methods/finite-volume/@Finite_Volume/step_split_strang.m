% Sat 11 Nov 12:21:35 CET 2017
% see ch 17.4 in randall and leveque
%
%% step in time, treat inhomogeneous part by Strang splitting 
%% this scheme is not suitable for stationary solutions, for example
%% steady shallow water flow
%
% TODO: the boundary condition in the advection step needs to be adapted,
%       to get a second order accurate solution
function [q obj] = step_split_strang(obj,t,q,dt)

	p = 0.5;

%	figure(1)
%	clf
%	plot(q)
%	pause

%	bcfun = @(t,q) obj.apply_bc(t,q);

	% 1/2 time step with source term
        q = obj.apply_bc(t, q);
	q = obj.odestepper(t, q, @(t,y) obj.sourcefun(t,obj.x,y,p*dt), p*dt);

%	figure(1)
%	clf
%	plot(q)
%	drawnow
%	pause(0.2)


	% full time step of advection (rieman problem)
        q = obj.apply_bc(t, q);
	q = obj.advect(t,q,dt);

	% 1/2 time step with source term
        q = obj.apply_bc(t+p*dt, q);
	q = obj.odestepper(t+p*dt, q, @(t,y) obj.sourcefun(t,obj.x,y,(1-p)*dt), (1-p)*dt);
end

%function [y obj] = step_split_strang(obj,t,x,y,f,g,dt)
%	% 1/2 time step with source term
%	y = rk2_step(t,x,y,f,0.5*dt);
%	% full time step of advection (rieman problem)
%	y = g(t,x,y,dt);
%	% 1/2 time step with source term
%	y = rk2_step(t,x,y,f,0.5*dt);
%end
