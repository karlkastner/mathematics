% 2015-01-07 14:38:42.552442646 +0100
% Karl Kastner, Berlin
%
%% make matrix A orhogonal to B
function [A, B] = orth2(A,B)
	A = orth(A);
	B = orth(B);
%	B'*A
	for idx=1:size(B,2)
		n = B(:,idx)'*A;
		A = A - B(:,idx)*n;
		% renormalise
		for jdx=1:size(A,2)
%			A(:,jdx)'*A(:,jdx)
			A(:,jdx) = 1/sqrt(A(:,jdx)'*A(:,jdx))*A(:,jdx);
%			A(:,jdx)'*A(:,jdx)
		end
	end
	A = orth(A);	
end

