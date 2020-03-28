% 2010-11-07 11:23
% Karl Kästner, KTH Stockholm

% TODO : port to java
% TODO, exploit symmetry (if matrix is symmetric)

% Bvec = [ a12 a23 ..  an-2n-1 an-1n          0;
%          a11 a22 a33      .. an-1n-1      ann;
%            0 a21 a32     a43      ..  a_n_n-1]
% A = [ [diag(B,+1).' 0]; diag(B).'; [0 diag(B,-1).']]
% syms a11 a12 a21 a22 a23 a33 a34 a35 a45                    
% B = [a11 a12 0 0; a21 a22 a23 0; 0 a33 a34 a35; 0 0 a44 a45]
%B = full(spdiags([ [9 10 11 12]' [5 6 7 8]' [1 2 3 4]'],-1:1,4,4))
%A = [ [diag(B,+1)' 0]; diag(B)'; [0 diag(B,-1)']]
%A = full(spdiags(A,-1:1,4,4))

function [B u_new] = qr_implicit_givens_2(B, varargin)
	if (nargin() > 1)
		shift = varargin{1};
	else
		shift = 0.0;
	end

	% bandwidth speed of for banded matrices
	if (nargin() > 2)
		error('','banded matrices are not supported');
	end
%		band = varargin{2};
%	else
%		band = 0;
%	end % band

	% quick check for hessian form
	if (length(B) > 2 && abs(B(3,1)) > 32*eps)
		throw('input matrix is not hessian');
	end
	% allocate the vector to store multipliers
	qV = zeros(size(B,1)*2,2);

	% matrix to vector
	A = [ [diag(B,+1); 0] diag(B) [0; diag(B,-1)] zeros(size(B,1),1)];

	% do the QR=A => A=RQ step
        [A u_new] = qr_implicit_givens_inner(A, qV, shift);

	% reconstruct the matrix
	B = spdiags(A,-1:2, size(B,1), size(B,1));
end % qr_implicit_givens_2

function [A u_new] = qr_implicit_givens_inner(A, qV, shift)
	n = size(A,1);

	% shift support
	if (0 ~= shift)
		u_old = shift;
		tmp = [A(n-1,2) A(n-1,1)
                       A( n,3)  A(n,2)];
		[uu u_new] = shiftQR(tmp, u_old);
		% apply shift B := B - uu*I;
 		% TODO, this can be done implicitely in the code  below
		A(:,2) = A(:,2) - uu;
	else
		uu = 0;
		u_new = 0;
	end % shift

	% first part : R = Q*B
	for idx=1:n-1
		alpha = A(idx,1);
		beta  = A(idx,2);
		k     = sqrt(abs(beta)^2 + alpha^2);
		u     = beta/k;
		v     = alpha/k;	
		q     = [u' v; -v u];
		qV(idx*2-1:idx*2,:) = q;
		% index for banded matrices
%		r = n*(~band) + min(n,(idx+2))*band;
%		B(idx:idx+1,idx:r) = q*B(idx:idx+1,idx:r);
% TODO band, for now just 1
		tmp = [A(  idx,2) A(idx+1,3) 0;
                        A( idx,1) A(idx+1,2) A(idx+1,1)];
		tmp = q*tmp;
		A(idx,2)   = tmp(1,1);
		A(idx+1,3) = tmp(1,2);
		A(idx+2,4) = tmp(1,3);
		A(idx , 1) = tmp(2,1);
		A(idx+1,2) = tmp(2,2);
		A(idx+2,3) = tmp(2,3);
	end % for idx

	% last row
	beta = A(n,2);
	k = abs(beta);
	u_end = beta/k;
%	B(end,end) = u_end'*B(end,end);
	A(n,2) = u_end'*A(n,2);

	% second part : R*Q' = Q*B*Q'
	for idx=1:n-1
		q = qV(idx*2-1:idx*2,:);
% TODO, band is 1 for now
%		t = 1*(~band) + max(1,(idx-1))*band;
%		B(t:idx+1,idx:idx+1)*q';
%		B(t:idx+1,idx:idx+1) = B(t:idx+1,idx:idx+1)*q';

		if (1 == idx)
		tmp = [A(  idx,2) A(idx+1,3)
                       A(  idx,1) A(idx+1,2)];
			tmp = tmp*q';
			A(  idx,2) = tmp(1,1);
			A(idx+1,3) = tmp(1,2);
			A(idx+1,2) = tmp(2,2);
			A(  idx,1) = tmp(2,1);

		else
		tmp = [A(  idx,3) A(idx+1,4);
                       A(  idx,2) A(idx+1,3);
                       A(  idx,1) A(idx+1,2)];
			tmp = tmp*q'
			A(  idx,3)   = tmp(1,1);
			A(idx+1,4)   = tmp(1,2);
			A(  idx,2)   = tmp(2,1);
			A(idx+1,3)   = tmp(2,2);
			A(idx+1,2)   = tmp(3,2);
			A(  idx,1)   = tmp(3,1);
		end
	end % for idx

	% last row
	A(n,3) = A(n,3)*u_end;
%	B(:,end) = B(:,end)*u_end;

	% undo shift
	A(:,2) = A(:,2) + uu;
%	B(1:n+1:n^2) = B(1:n+1:n^2) + uu;
	% do not force B to be symmetric (effort 0.5*n^2 !!!)
%	A - full(B)
%	full(spdiags(A,-1:2,n,n)) 
%		pause
%	full(A)
end % qr_implicit

