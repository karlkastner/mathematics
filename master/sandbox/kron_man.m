function C = kron_man(A, B)
	ma = size(A,1);
	na = size(A,2);
	mb = size(B,1);
	nb = size(B,2);
	C = zeros(ma*mb,na*nb);	
	for idx=1:size(A,1)
		for jdx=1:size(A,2)
			C( (idx-1)*mb+1:idx*mb, (jdx-1)*nb+1:jdx*nb ) = A(idx,jdx)*B;
		end
	end
end

