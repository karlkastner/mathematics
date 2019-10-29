function C = ichol_man(A)
	% allocate memory
	C = triu(A);
	for k=1:n
		C(k,k) = sqrt(C(k,k));
		id = 1.0/C(k,k);
		C(k,k+1:end) = id*C(k,k+1:end);
		for j=find(C(k+1:n))); %k+1:n
		 for i=j:n
			
			A(i,j) = A(i,j) - A(i,k)*A(j,k);
		 end
		end
	end
end

