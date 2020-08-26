% Sat 28 Oct 14:43:06 CEST 2017
function [AA,bb] = assemble_AA(obj,cdx,ypm)

	nxc    = obj.nxc(cdx);
	nci    = obj.nci(:,cdx);
	npi    = obj.npi(:,cdx);

	% construct function values at section mid-points from separated parts
	yc = zeros(nci(end)-1,1);
	for edx=1:obj.neq
		switch (obj.oo(edx))
		case {1}
			% first order ode
			%c = cc(1,end-2:end,edx);
			% TODO this condition has to be checked for each row (!)
			if (false) %0 ~= c(1,2)) 
			yc(nci(edx):nci(edx+1)-1) = ( ...
			   		  ypm(npi(edx)  :2:npi(edx+1)-2) ...
			   		+ ypm(npi(edx)+1:2:npi(edx+1)-1) ...
					);
			else
			if (~obj.opt.dischargeisvariable)
				yc(nci(edx):nci(edx+1)-1) = ... 
				   		ypm(npi(edx)+1:2:npi(edx+1)-1);
			else
				yc(nci(edx):nci(edx+1)-2) = ... 
				   		ypm(npi(edx)+1:2:npi(edx+1)-2);
				yc(nci(edx+1)-1) = ypm(npi(edx+1)-1);
			end	
			end

		case {2}
			% 2nd order ode
			% sum of left and right going wave as well as inhomogeneous part
			yc(nci(edx):nci(edx+1)-1) = ( ...
					  ypm(npi(edx)  :3:npi(edx+1)-3) ...
					+ ypm(npi(edx)+1:3:npi(edx+1)-2) ...
				 	+ ypm(npi(edx)+2:3:npi(edx+1)-1) ...
				 	);
		otherwise
			error('here');
		end	
	end % for edx

	% ode-coefficients at section mid-points
	 obj.out(cdx).cc = feval(obj.odefun,cdx,obj.out(cdx).xc,yc);

	% TODO, global buffer
	AA = sparse(npi(end)-1,npi(end)-1,6*npi(end));
	bb = zeros(npi(end)-1,1);
	obj.out(cdx).ll = zeros(nxc,2,obj.neq);

	% for each of the coupled odes
	for edx=1:obj.neq
		% eigenvalues, roots of the characteristic polynomial
		% roots are in general not complex conjugate pairs,
		% as this is only the case when coefficients are real
		% and case the solution a damped wave
		obj.out(cdx).ll(:,:,edx) = roots2(obj.out(cdx).cc(:,1:3,edx));
%		ll(:,:,edx) = l;

		% asumption : ode are either purely 1st or 2nd order, but not mixed
		switch (obj.oo(edx))
		case {1}
			if (obj.opt.dischargeisvariable)
				[A, b] = obj.assemble1_A_Q(cdx,edx);
			else
				[A, b] = obj.assemble1_A(cdx,edx);
			end
		case {2}
			[A, b] = obj.assemble2_A(cdx,edx);
		otherwise
			error('here');
		end
	
		% stack system of odes
		% when the equations are only weekly non-linear,
		% so they can be solved individually
		AA(npi(edx):npi(edx+1)-1, npi(edx):npi(edx+1)-1) = A;
		bb(npi(edx):npi(edx+1)-1) = b;
	end % for edx
%	obj.out(cdx).ll = ll;
%	obj.out(cdx).cc = cc;
end


