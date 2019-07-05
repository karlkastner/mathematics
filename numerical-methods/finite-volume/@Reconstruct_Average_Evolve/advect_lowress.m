% 2011-02-20 15:55:20.000000000 +0100
%% single time step
%% low resolution
function q = step(obj, t, q, dt)
	qold = q;	
if (0)
% variant 1
        f = feval(obj.flux, t, q);
        q_ = q - dt/obj.dx*obj.Dr*f;
        q_ = feval(obj.bc, t, q_); % necessary ?
        f_ = feval(obj.flux, t+dt/2, q_); %, g, M);
        q__ = q_ - dt/obj.dx*obj.Dl*f_;
        qv1 = 0.5*(q + q__);
        % apply boundary condition
        qv1 = feval(obj.bc, t, qv1);
else
	% variant 2 - almost ident
	% Richtmeyr two step lax wendroff
	% TODO alternative implementation had 1 dt and 1/2 dt
	% in half steps instead of 1/2 dt and 1 dt

	% compute fluxes at volume interfaces inbetween grid points
        f = feval(obj.flux, t, q);
	% Randall Leveque eq. 4.23
	% half-step with left and right differences	
	% note: 0.5 in m was originally missing
	ql = 0.5*obj.Ml*q - 0.5*dt/obj.dx*obj.Dl*f;
	qr = 0.5*obj.Mr*q - 0.5*dt/obj.dx*obj.Dr*f;

	force = zeros(2*length(obj.x),1);
	for idx=1:length(obj.forcefun)
		%force = force + obj.forcefun(t, obj.x, q);
		force = force + obj.forcefun{idx}(t, obj.x, q);
	end
	% TODO left and right differences in the force?
	%norm(force)
	ql = ql + 0.5*dt*force;
	qr = qr + 0.5*dt*force;
	
	% compute fluxes at volume midpoints at the grid points
        Fl = feval(obj.flux, t+0.5*dt, ql);
        Fr = feval(obj.flux, t+0.5*dt, qr);

	% Randall-Leveque Eq. 4.4
	% half step with central differences
	% FTCS
	q = q - dt/obj.dx*(Fr - Fl);

	% TODO full step here?
	q = q + dt*force;		

        % apply boundary condition
        q = obj.apply_bc(t, q);

	%max(q(1:end/2))
	%dt	
	%figure(1)
	%plot(q(1:end/2))
	%dt
	%pause(0.1)
end

