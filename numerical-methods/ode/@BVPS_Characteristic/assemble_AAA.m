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
%function [AAA,bbb] = assemble_AAA(obj,ypm)
%
function [AAA,bbb] = assemble_AAA(obj,ypm)
	npii = obj.npii;

	% allocate memory for differential operator
	% TODO global buffers
	AAA  = sparse([],[],[],npii(end,end)-1,npii(end,end)-1,6*npii(end,end)-1);

	% allocate memory for inhomogeneous part
	bbb  = zeros(npii(end,end)-1,1);

	% for each edge in graph (channel in network)
	for cdx=1:obj.nc
		% set up discretisation matrix for coupled odes along edge (channel)
		[AA,bb] = obj.assemble_AA( cdx, ypm(npii(1,cdx):npii(end,cdx)-1) );

		% write to global discretization matrix
		AAA(npii(1,cdx):npii(end,cdx)-1,npii(1,cdx):npii(end,cdx)-1) = AA;
		bbb(npii(1,cdx):npii(end,cdx)-1,1)     = bb;
	end % for cdx

	if (~isempty(obj.jfun))
		[AAA,bb] = obj.couple_junctions(AAA,bb);
	end

end % assemble_AAA

