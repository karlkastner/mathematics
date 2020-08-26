% Sat 28 Oct 14:43:06 CEST 2017
% Karl Kastner, Berlin
%
%% solve system of non-linear second order odes (in more than one variable)
%% as boundary value problems
%%
%% odefun provides ode coefficients c:
%% c(x,1) y''(x) + c(x,2) y'(x) + c(x,3) y = c(x,4)
%%    c_1 y"     + c_2 y'       + c_3 y + c_4 = c_4
%%
%% subject to the boundary conditions
%% bcfun provides v and p and optionally q, so that:
%%
%% b_1 y + b_2 y' = f
%%    q(x,1)*( p(x,1) y_l(x) + p(x,2)  y_l'(x)
%%  + q(x,2)*( p(x,1) y_r(x) + p(x,2) y_r'(x)    = v(x)
%% where q weighs the waves travelling from left to right and right to left (default [1 1])
%
% TODO allow for resolve (do not rerun through init process)
function solve(obj)

	ni = obj.ni;
	nci = obj.nci;
	npi = obj.npi;
	npii = obj.npii;

	% TODO move to (re)-init ?
	% initial value of ypm
	% complex amplitude of the left and right going wave at segment mid points
	ypm    = zeros(obj.npii(end,end)-1,1);
	if (~isempty(obj.inifun))
	    for cdx=1:obj.nc
		yi      = feval(obj.inifun,cdx,obj.out(cdx).x);
		if (length(yi) == npi(end,cdx)-1)
				ypm(npii(1,cdx):npii(end,cdx)-1) = yi;
		else
		    for edx = 1:obj.neq
			% assign initial values into the inhomogeneous part
			% this does not matter, as yi is reassembled later
			% the inhomogeneous part is always the second part
			if (length(yi) == nci(end,cdx)-1)
				yi_ = yi(nci(edx,cdx):nci(edx+1,cdx)-1);
			else
				yi_ = yi(ni(edx,cdx):ni(edx+1,cdx)-1);
			end
			switch (obj.oo(edx))
			case {1}
				if (~obj.opt.dischargeisvariable)
					if (length(yi) ~= nci(end,cdx)-1)
						yi_ = mid(yi_);
					end
					ypm(npii(edx,cdx)+1:2:npii(edx+1,cdx)-1) = yi_;
				else
					% discharge
					ypm(npii(edx+1,cdx)-1) = yi_(end);

					% level
					if (length(yi) ~= nci(end,cdx)-1)
						yi_ = mid(yi_(1:end-1));
					end
					ypm(npii(edx,cdx)+1:2:npii(edx+1,cdx)-2) = yi_;
				end
			case {2}
				if (length(yi) ~= nci(end,cdx)-1)
					yi_ = mid(yi_);
				end
				ypm(npii(edx,cdx)+1:3:npii(edx+1,cdx)-1) = yi_;
			end % switch
		   end % if
		end % for edx
	    end % for cdx
	end % ~isempty(ifun)

	% solve non-linear system by picard iteration
	% TODO iteration results should not be stored in out(1)
	[ypm, obj.out(1).cflag, obj.out(1).kiter] = picard(@solve_,ypm,obj.opt.sopt);

	% interpolate solution to segment end points (grid points)
	obj.reconstruct(ypm);

function ypm = solve_(ypm)
	[A,b] = obj.assemble_AAA(ypm);

	% balance
	if (obj.opt.balance)
		% todo, two sided with sqrt(s)
		s  = 1./sqrt(abs(diag(A)));
		A = diag(s)*A.*diag(s);
		b = s.*b;
	end

	% solve
	ypm    = (A \ b);
%	ypm = pinv(full(AA))*bb;
end % solve_

end % function solve

