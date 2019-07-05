% Sa 6. Feb 23:38:49 CET 2016
% Karl Kastner, Berlin
%% reconstruct matrix from eigenvalue decomposition
function A = ieig(V,E)
	E  = diag2x2(E);
	VE = mtimes2x2(V,E);
%	VE = zeros(size(V));
%	VE(1,1,:) = V(1,1,:).*E(1,:);
%	VE(1,2,:) = V(1,2,:).*E(2,:);
%	VE(2,1,:) = V(2,1,:).*E(1,:);
%	VE(2,2,:) = V(2,2,:).*E(2,:);

	Vi = inv2x2(V);

	A = mtimes2x2(VE,Vi);
%	A = zeros(size(V));
	% V E V'
%	A(1,1,:) = VE(1,1,:).*Vi(1,1,:) + VE(1,2,:).*Vi(2,1,:);
%	A(1,2,:) = VE(1,2,:).*Vi(1,2,:) + VE(1,2,:).*Vi(2,2,:);
%	A(2,1,:) = VE(2,1,:).*Vi(1,1,:) + VE(1,2,:).*Vi(2,1,:);
%	A(2,2,:) = VE(2,1,:).*Vi(1,2,:) + VE(2,2,:).*Vi(2,2,:);
end

