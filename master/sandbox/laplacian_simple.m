% 2012 Dec  4 04:12 MSK
% Karl Kästner, Berlin

% laplacian with dirichlet boundary conditions
function L = laplacian_simple(n,order,mode)
	if (nargin < 2)
		order = 2;
	end
	if (nargin < 3)
		mode = 'strong';
	end
	switch (order)
		case {2}
			% 2nd order
			kernel = [1    -2     1];
			k_den = 1;
		case {4}
			% 4th order
			kernel = [-1   16  -30   16  -1 ];
			k_den = 12;
			ip     = [ 5  -10   10   -5   1 ];
		case {6}
			% 6th order
			kernel = [ 2   -27   270  -490   270   -27    2 ];
			k_den = 180;
			ip     = [ 7   -21    35   -35    21    -7    1 ;
			          28  -112   210  -224   140   -48    7 ];

		case {8}
			% 8th order
			kernel	= [  -9	 128 -1008  8064 -14350  8064 -1008  128  -9 ];
			k_den = 5040;
			ip      = [   9  -36    84  -126    126	  -84    36   -9   1 ;
			             45 -240   630 -1008   1050  -720   315  -80   9 ;
			            165 -990  2772 -4620   4950 -3465  1540 -396  45 ];
		case {10}
			% 10th order
			kernel	= [    8   -125   1000	-6000  42000  -73766  42000  -6000   1000   -125      8];
			k_den = 25200;
			ip      = [   11    -55    165   -330    462    -462    330   -165     55    -11      1 ;
			              66   -440   1485  -3168   4620   -4752   3465  -1760    594   -120     11 ;
			             286  -2145   7722 -17160  25740  -27027  20020 -10296   3510   -715     66 ;
			            1001  -8008  30030 -68640 105105 -112112  84084 -43680  15015  -3080    286 ];
		otherwise
			disp(o)
			'error'
			return;
	end
	% stepwidth
	h = 1/(n+1);
	% laplacian
	L = spdiags(ones(n,1)*kernel, -order/2:order/2, n, n);

	switch (mode)
	case {'central'}
		% substract from central (preserves symmetrie)
		%full(L(1:10,1:10))
		%pause
		%L
		for idx=1:order/2-1
			for jdx=1:order/2-idx % first is zero (homogenous dirichlet)
				L(idx,idx) = L(idx,idx) - kernel((jdx)); %+ abs(kernel(jdx));
		%full(L(1:10,1:10))
		%pause
			end
			L(end+1-idx,end+1-idx) = L(idx,idx);
		end
		%full(L(1:10,1:10))
		%pause
	case {'lower'}
		% lower order fdm at the boundary
		% approximate boundaries with lower order FDM (not as good as interpolation)
% h^2 wrong - applied below
% den wrong - undo
		if (order >= 4)
			L(1,:) = 0;
			L(1,1:2) = [-2 1]*k_den;%/h^2;
			L(end,:) = 0;
			L(end,end-1:end) = [1 -2]*k_den;%/h^2;
		end
		if (order >= 6)
			L(2,:) = 0;
			L(2,1:4) = [16  -30   16  -1 ]*k_den/12*1;%/h^2;
			L(end-1,:) = 0;
			L(end-1,end-3:end) =  [-1   16  -30   16]*k_den/12*1;%/h^2;
		end
		if (order >= 8)
			L(3,:) = 0;
			L(3,1:6) = [-27   270  -490   270   -27    2 ]*k_den/180*1;%/h^2;
			L(end-2,:) = 0;
			L(end-2,end-5:end) = [2 -27   270  -490   270   -27]*k_den/180*1;%/h^2;
		end
		if (order >= 10)
			L(4,:) = 0;
			L(4,1:8) = [128 -1008  8064 -14350  8064 -1008  128  -9 ]*k_den/5040 * 1;%/h^2;
			L(end-3,:) = 0;
			L(end-3,end-7:end) = [  -9  128 -1008  8064 -14350  8064 -1008  128]*k_den/5040 * 1;%/h^2;
		end
	case {'weak'}
		% weakly impose the boundary condtions
		%L(1,:) = 0;
		%L(end,:) = 0;
%		for idx=1:order/2
%			L(idx,idx) = abs(L(idx,idx))*1e7;
%			L(end-idx+1,end-idx+1) = abs(L(end-idx+1,end-idx+1))*1e7;
%		end
		L(order/2,order/2) = abs(L(order/2,order/2))*1e7;
		L(end+1-order/2) = abs(L(end+1-order/2,end+1-order/2))*1e7;
		%h = 1/(n-order/2);
		h = 1/(n+1-order);
	case {'strong'}
		% extrapolate boundary points, dirichlet boundary, e.g. first outer point is zero
		% todo implement as matrix-matrix multiplication
		for rdx=1:order/2-1
		 for cdx=1:order/2-rdx
			L(rdx,1:order)               = L(rdx,1:order)               + kernel(cdx)*ip(1-rdx+order/2-cdx,2:end);
			L(end-rdx+1,end-order+1:end) = L(end-rdx+1,end-order+1:end) + kernel(cdx)*fliplr(ip(1-rdx+order/2-cdx,2:end));
		 end
		end
	case {'none'}
		% nothing to do
	otherwise
		'error'
	end % switch

	% apply stepwidth
	L = 1/k_den*1/h^2*L;
end % laplacian_simple

