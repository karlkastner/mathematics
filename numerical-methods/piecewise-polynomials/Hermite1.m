% Sun May  1 15:50:24 CEST 2016
% Karl Kastner, Berlin
%% hermite polynomial interpolation in 1d
classdef Hermite1

	methods (Static)

	% parametrisation from -1,0,1, fourth parameter given by constrained derivative
	function c = set_up(f,ci)
%		f = (1:9)';
		f = rvec(f);

		n = length(f);
		nb = ceil(length(f)/(3-ci));
		switch (ci)
		case {0}
			f(end+1:end+2) = NaN;
			n = n+2;	
		case {1}
			% padd
			f(end+1) = NaN;
%			n = n+1;
%			r  = mod(n,(3-ci));
%			f(end+1:end+r) = f(end);
%			n = n+r;
		case {2}
		end

		switch (ci)
		case {0}
			H = 	[1  -1  1  -1;
			         1   0  0   0;
			         1   1  1   1;
				 1   2  4   8];  
	
			% set up block matrix for piecewise polynomials
			I  = speye(nb);
			HH = kron(I,H);
			A  = HH;
	
			% set up rhs
			rhs = [f(1:3:n-3);f(2:3:n-2);f(3:3:n-1);f(4:3:n)];
			rhs = rhs(:);
				
			%A = A(1:n-1,1:n-1);
			%rhs = rhs(1:end-1);
			rhs = rhs(1:n+nb-3);
		case {1} % hermite	
% -> continuity for non-perfect block split is problematic,
%    make partial block a full block by overlapp with previous block	
			% hermite block
			H = 	[0, -1, 2, -3;
			         1  -1  1  -1;
			         1   0  0   0;
			         1   1  1   1];
			C = [ zeros(3,4);
			      0 1 2 3 ];
	
			% set up block matrix for piecewise polynomials
			I  = speye(nb);
			HH = kron(I,H);
			CC = kron(I,C);
			CC = [zeros(1,(nb)*4); CC(1:end-1,:)];
			A = HH + CC;
	
			% set up rhs
			rhs = [f(1:2:2*(nb+1)-2);f(2:2:2*(nb+1)-1);f(3:2:2*(nb+1)); zeros(1,nb)];
			rhs = rhs(:);

			A(1,:) = [];
			%rhs = rhs(1:end-1);	
			%rhs = rhs(1:n+2*nb-2);
			rhs = rhs(1:n+2*nb-2);
		case {2}
			nb = n;
			H = 	[0, -1,  2, -3;
				 0,  0, -2,  6;
			         1  -1  1  -1;
			         1   0  0   0];
			C = [ zeros(2,4);
			      0 1 0 0
	                      0 0 2 0 ];
	
			% set up block matrix for piecewise polynomials
			I  = speye(nb-1);
			HH = kron(I,H);
			CC = kron(I,C);
			CC = [zeros(2,(nb-1)*4); CC(1:end-2,:)];
			A = HH + CC;
	
			% set up rhs
			rhs = [f(1:1:n-1); f(2:1:n); zeros(2,nb-1)];
			rhs = rhs(:);
	
			rhs = rhs(1:end-2);
			A(1:2,:) = [];
			A = A(:,1:end-2);
		otherwise
			error('here');
		end
			% delete excess rows and cols
			nr = length(rhs);
			A = A(1:nr,1:nr); 
		[full(A) rhs A\rhs]
ci
pause
%		size(A)
%		size(rhs)
		% determine coefficients
		c = A \ rhs;
		% stack into blocks
		r = mod(length(c),4);
		if (r>0) r = 4-r; end
		c(end+1:end+r) = 0;
		r = mod(length(c),4);
		c = reshape(c,4,[]);
	end


	% interpolation
	function fi = interpolate(xi)
		% determine block interval of target point
		block = bsxfun(@greater,x,x0);
		block = sum(block,2);
		block = min(nblock,max(1,block));
		% transform target coordinates to block interval coordinates
		x = 2*(x-x0(block))./l(block) - 1;
	
		for idx=1:np
			% set up the cubic vandermonde matrix
			A = vander_1d(x,3);
			% evaluate source polynomial
		end
	end

	function d = derivative1(f,ci);
		c   = Hermite1.set_up(f,ci);
		d   = Hermite1.derivative1_(c,ci);
			% make a vector
			d = reshape(d,[],1);
		d(end+1:end+2) = NaN;
		d = d(1:length(f));
	end

	function d = derivative2(f,ci);
		c   = Hermite1.set_up(f,ci);
		d   = Hermite1.derivative2_(c,ci);
			% make a vector
			d = reshape(d,[],1);
		d(end+1) = NaN;
		d = d(1:length(f));
	end

	function d = derivative1_(c,ci)
		% derivative matrix
			switch (ci)
			case {0}
			% derivative matrix
			D1 = 	[0,  1,-2, 3;
				 0   1  0   0;
			         0   1  2   3];
			case {1}
				% derivative matrix
				D1 = 	[0,  1,-2, 3;
					 0   1  0   0;];
			case {2}
				% derivative matrix
				D1 = 	[0,  1,-2, 3];
			otherwise
				error('here');
			end
			d  = D1*c;
	end
	function d = derivative2_(c,ci)
			% second derivative matrix
			switch (ci)
			case {0}
			D2 =   [0, 0, 2, -6;
				0, 0, 2,  0;
			        0, 0, 2,  3];
			case {1}
				D2 =   [0, 0, 2, -6;
					0, 0, 2,  0;];
			case {2}
				D2 =   [0, 0, 2, -6];	
			otherwise
				error('here');
			end
		d  = D2*c;
	end

	function curv = curvature(x,y,ci);
		cx   = Hermite1.set_up(x,ci);
		cy   = Hermite1.set_up(y,ci);
		curv = Hermite1.curvature_(cx,cy,ci);
			curv = curv(1:length(x));
	end
	% TODO for c0 and c1 (hermite) interpolation, the second derivative jumps between blocks, average left and right
		function curv = curvature_(cx,cy,ci)
			% at each point get dx/dt, dy/dt,d^x/dt^2 and d^y/dt^2
	
			% evaluate
			dx  = Hermite1.derivative1_(cx,ci);
			dy  = Hermite1.derivative1_(cy,ci);
			ddx  = Hermite1.derivative2_(cx,ci);
			ddy  = Hermite1.derivative2_(cy,ci);
		
			curv = (dx.*ddy - dy.*ddx)./(dx.^2 + dy.^2).^1.5;
			
			% make a vector
			curv = reshape(curv,[],1);
		end % curvature_
	end % methods (Static)
end % class Hermite1

