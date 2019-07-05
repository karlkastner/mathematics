% Thu Apr 28 04:51:19 MSD 2011
% Karl Kästner

%% forward euler method with staggered grid
function u = stagger_euler(t, u, bc)
        Dh = spdiags(ones(n-1,1)*[-1 1], 0:1, n-1, n);
        Dv = spdiags(ones(n,1)*[-1 1], -1:0, n, n-1);
        Dvs = spdiags(ones(n-1,1)*[-0.5 0 0.5], -1:1, n-1, n-1);
        Dvs(1,:) = 0; Dvs(end,:) = 0;

        % h staggered
        h_ = u(1:n);
        v_ = u(n+1:end-1);
        hs = 0.5*(h_(1:end-1) + h_(2:end));
        
        h = h_ - 0.5*dt/dx*(Dv*(hs.*v_));
        v = v_ - 0.5*dt/dx*(0.5*(Dvs*(v_.*v_)) + g*Dh*h_);
% linearised
%        h = h_ - 0.5*dt/dx*(Dv*(v_));
%        v = v_ - 0.5*dt/dx*(g*Dh*h_);
        u = [h; v; 0];

        % apply bc
        u(1:end-1) = feval(bc, t, u(1,end-1));
end % staggered_euler

