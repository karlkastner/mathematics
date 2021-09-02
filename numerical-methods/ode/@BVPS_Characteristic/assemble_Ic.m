% Mon 14 Sep 11:36:32 +08 2020

function assembe_Ic(obj)
	% TODO allocate
	Abuf = [];
	for cdx=1:obj.nc
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
			id = nbuf+(1:(nci(edx+1)-nci(edx)));
%			Abuf(,1) = nci(edx):nci(edx+1)-1;
%			Abuf(,2) = (npi(edx)  :3:npi(edx+1)-3);
%			Abuf(,3) =  1
%			nbuf = nbuf+
%			Abuf(,1) = nci(edx):nci(edx+1)-1;
%			Abuf(,2) = (npi(edx)+1  :3:npi(edx+1)-2);
%			Abuf(,3) =  1
%			nbuf = nbuf+
%			Abuf(,1) = nci(edx)+1:nci(edx+1)-2;
%			Abuf(,2) = (npi(edx)+2:3:npi(edx+1)-1);
%			Abuf(,3) =  1
			%nbuf = nbuf+
		otherwise
			error('here');
		end	
	end % for edx
	end % for cdx

	obj.Ic = sparse(Abuf(:,1),Abuf(:,2),Abuf(:,3),nxc_,nxc_);
end

