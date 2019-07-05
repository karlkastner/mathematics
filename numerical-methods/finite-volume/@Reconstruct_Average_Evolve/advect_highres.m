% 2011-02-20 15:55:20.000000000 +0100
% Sun 12 Nov 15:20:45 CET 2017
% Karl Kastner, Berlin
%% single time step for the reconstruct evolve algorithm
% TODO homogenization in case of rapidly varying coefficients may be necessary (cf. ch 9.14)
% TODO entropy fix for 0==s for hydraulic jumps
% TODO, use lax-liu method instead of 9.75
function [q obj] = step(obj,t,q,dt)
	dx   = obj.dx;
	m    = obj.pde.m;
	n    = length(obj.x);

	% matrix to vector
	m2w  = kron(eye(m),ones(n,1));

        % apply boundary condition
        q = obj.apply_bc(t, q, dt);

	% unlimited LW-update:
	% u = u - 0.5 dt/dx (J ur - J ul) 

	% 1.) roe-average at left interface
	%     fr = J Dr q
	[qhat auxhat] = obj.pde.roe_average(q,obj.Left*q);

	% 2.) eigendecomposition of jacobian flux matrix
	% TODO jacobian or true matrix?
	%     fr_ = iV fr = L iV Dr q
	% randalls notation:  Vl L Vr D q = sum l_p vl_p vr_p dq = sum l_p w_p, where s = l in the nonlinear case
	[L Vr Vl] = obj.pde.fluxmateig(qhat,auxhat);
	l  = diag(L);
	l1 = l(1:n);
	l2 = l(n+1:2*n);
	s  = [[l1; l1],[l2;l2]];

	% differences of q (not qhat)
	Dlq = obj.Dl*q;

	% upwind diffences of q for all eigenvalues, 6.61
	% wave, 9.72
	% w = a vr = vl dq vr (6.55), 9.75, 3.23
	w_l = Vr*bsxfun(@times,Vl*Dlq,m2w);

	w_up = (s>0).*(obj.Left*w_l) + (s<0).*(obj.Right*w_l);

	% TODO entropy fix

	% compute input value for flux limiter, 9.74
	% since the waves are not co-linear, comparing scalar scales of the vectors is not consistent
	% for linear case theta = aa_up./aa;
	theta = dot(w_up,w_l)./dot(w_l,w_l);
	theta(~isfinite(theta))=0;
	
	% limit wave, 9.69
	w_l_tilde = repmat(obj.limiter(theta),m,1).*w_l;

	% 3.) apply flux limiter, so that
	%     fr_ = (1-phi) fr_up + phi fr_lw
	% 5.) recompute fluxes
	%     fr = V fr

	% 9.67 limited flux through left interface
	% for linear schemes this simplifies to (c.f 6.9):
	% Fl := 1/2 (|A| - dt/dx*A^2) = Vl*(lambda-dt/dx*lambda^2)*Vr
	Fl = 0.5*sum(abs(s).*(1-dt/dx*abs(s)).*w_l_tilde,2);
	% flux through right interface
	%Fr = obj.Right*Fl;

	AmD = obj.Right*sum((s.*(s<0)).*w_l,2);
	ApD = sum((s.*(s>0)).*w_l,2);

	% 6.) time step update 9.66
	% NB: Am Dlq + Ap Dlq in 6.59 vs Ap Drq + Am Dlq in 9.66
	% A Dq = s W = lambda W
	%q = q + dt/dx*(Ap*(Drq) + Am*(Dlq)) ...
	%frtilde = (1-phi(r))fupr + phi(r)flwr

	Fr = obj.Right*Fl;

	% source term (reactive only)
%	xc     = leftmean(obj.x);
	f_l    = obj.pde.sourcefun(t,obj.x,qhat,auxhat);

%	g      = 9.81;
%	zb     = obj.zbfun(obj.x);
%	Dlzb   = ldiff(zb)/dx;
%	h_hat  = qhat(1:n);
%	f_l    = [zeros(n,1);-g*h_hat.*Dlzb];
%	norm(f_l-f_l_)
%pause

	% transform
	f_l    = Vr*bsxfun(@times,Vl*f_l,m2w);
        fm     = obj.Right*sum(((s<0)).*f_l,2);
        fp     =           sum(((s>0)).*f_l,2); 

	% not sure if this is right
	f_l_tilde = repmat(obj.limiter(theta),m,1).*f_l;
	Sl        = 0.5*sum(abs(s).*(1-dt/dx*abs(s)).*f_l_tilde,2);
	Sr        = obj.Right*Sl;

%	k=3;
%	Fl([1:k, n-k+1:n,n+1:n+k,2*n-k+1:2*n]) = 0;
%	Fr([1:k, n-k+1:n,n+1:n+k,2*n-k+1:2*n]) = 0;
	
	q = q - dt/dx*(ApD + AmD) ...
	      - dt/dx*(Fr-Fl) ...
              + dt*(fm+fp);
	      + dt*(Sr-Sl);
		;
	      %- dt/dx*(obj.Dr*Fl);
	
	function ab = dot(a,b)
		%ab = a;
		for idx=1:size(a,2)
			ab(:,idx) = sum(reshape(a(:,idx).*b(:,idx),[],m),2);
		end
	end
end % step

