% Wed Jul 11 19:46:44 MSK 2012
% Karl Kästner, Berlin

function C_dx = derivarive_1d(C,n)
	C_dx = zeros(n,size(C,2));
	for idx=1:n
		C_dx(idx,:) = idx*C(idx+1,:);
	end
end % derivative_1d

