% 2018-05-25 16:29:55.529134236 +0200 linear-algebra/orth_man.m
%% orthogonalize the columns of A
function A = orth_man(A)
	for idx=1:size(A,2)
		A(:,idx) = A(:,idx)/norm(A(:,idx));
	end
	for idx=1:size(A,2)-1
		for jdx=idx+1:size(A,2)
			A(:,jdx) = A(:,jdx) - (A(:,idx)'*A(:,jdx))*A(:,idx);
		end
	end
end
