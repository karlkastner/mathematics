function A = orth2(A,B)
	for idx=1:size(B,2)
		for jdx=1:size(A,2)
			A(:,jdx) = A(:,jdx) - (A(:,jdx)'*B(:,idx))*A(:,idx);
		end
	end
end
