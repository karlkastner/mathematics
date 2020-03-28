% 2011-12-19 09:33:35 UTC
% Karl Kästner, Berlin

% this seems to be an outdated duplicate of derive_fdm

% one sided approximation in the corner

function [dF dF_num dF_den] = derive_bc_one_sided()
n=3;
%switch(n)
%	case {1}
		A = -[-1  1  0;
      		       0  1 -1];
		B =  [-1  1 -1;
		       1  1  1];

		A = -[-1  1  0;
                       0  1  0;
      		       0  1 -1]
		B =  [-1  1 -1
                       0  0  0
                       1  1  1]
%{
		A = -[ 1 -1 0;
		       1  0 -1];
		% columns are powers of dx (dx = [1h 2h 3h ...]')
		dx = (1:n-1)'
		B = repmat(dx,1,n).^repmat(1:n,n-1,1)
		% coefficients of the Taylor expansion
		% c = 1/n!
%}
		C = diag(1./factorial(1:n))
%	case {2}
	A__  = [
	    -1     1     0     0     0
	     0     1     0     0     0
	     0     1    -1     0     0
	     0     1     0    -1     0
	     0     1     0     0    -1 ];
	% vandermonde matrix, dx^n
	% x = [-1 : 3];
	B__  = [
	     -1     1    -1     1
	      0     0     0     0
	      1     1     1     1
	      2     4     8    16
	      3     9    27    81 ];

A_= [
    -1     0     0     1     0
     0    -1     0     1     0
     0     0    -1     1     0
     0     0     0     1     0
     0     0     0     1    -1 ];
B_= [
     -3     9    -27    81
     -2     4     -8    16
     -1     1     -1     1
      0     0      0     0
      1     1      1     1 ];
A___ = [
    -1     0     1     0     0
     0    -1     1     0     0
     0     0     1     0     0
     0     0     1    -1     0
     0     0     1     0    -1 ];
B___ = [-2     4    -8    16
     -1     1    -1     1
      0     0     0     0
      1     1     1     1
      2     4     8    16
];


% A F = B C dF
% F = [f(-h) f(0) f(+h) f(2h) f(3h)]
B*C
dF = (B*C) \ A
dF_den = 1./min(abs(dF'))'
dF_num = dF .* (dF_den*ones(1,size(A,2)))

end % function derive_bc_one_sided

