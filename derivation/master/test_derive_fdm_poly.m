% 2012-04-19 08:18:40
% Karl Käsnter, Berlin

function test_derive_fdm_poly()

syms xk xl xr xc xr xs h dx
syms fk fl fr fc fr fs

if (0)
	% constant stepwidth
	X_ = (-2:2)*dx;
else
	X_ = [xk xl xc xr xs];
end

F_ = [fk fl fc fr fs]; 
%X_ = (-2:2)*h;
for jdx=1:2;
	X = X_(3-jdx:3+jdx);
	F = F_(3-jdx:3+jdx);
	for xc_=-jdx:0
		xc=X(xc_+3);
		syms T;
		T(1,2:1+length(X)) = F;
		for idx=2:length(X)
			T(idx,1) = idx-1;
			D  = derive_fdm_poly(X,xc);
			D_ = D(idx,:);
			%[num, den]=numden(D_*F.');
			%	T(idx-1,2) = simple(factor(den));
			%	T(idx-1,3:7) = simple(factor(D_*den)); %, (factor(den));
			%T(idx,2:1+length(X)) = factor(D_);
			T(idx,2:1+length(X)) = (D_);
		end % idx
		fprintf('xc = %d\n',xc_);
		T
	end % x_c
end % jdx

end % func test_derive_fdm_poly()

