% Sat 28 Oct 14:43:06 CEST 2017
%
% function [AAA,bbb,out] = bvp2c_assemble(out,ypm,odefun,bcfun,xi,neq,nxc,oo,ni,nci,npi,npii,jfun,dischargeisvariable)
%
%
% output :
%	AAA : left hand side discretization matrix
%	bbb : right hand side vector
%	out(id).lll : eingevalues of the characteristic functions
%
% function assemble_AAA(obj,ypm)
%
function [A,b] = assemble_AAA(obj,ypm)
	npii = obj.npii;


	% allocate memory for inhomogeneous part
	obj.b  = zeros(npii(end,end)-1,1);
	% TODO allocate
	obj.Abuf = [];
	obj.nbuf = 0;

	% for each edge in graph (channel in network)
	for cdx=1:obj.nc
		% set up discretisation matrix for coupled odes along edge (channel)
		obj.assemble_AA( cdx, ypm(npii(1,cdx):npii(end,cdx)-1) );

	end % for cdx

	A = sparse(obj.Abuf(:,1),obj.Abuf(:,2), obj.Abuf(:,3), obj.npii(end,end)-1, obj.npii(end,end)-1);
	b = obj.b;

	if (~isempty(obj.jfun))
		[A,b] = obj.couple_junctions(A,b);
	end

end % assemble_AAA

