% Sat Jan 21 02:35:26 MSK 2012
% Karl Kästner

function [coeff residual] = derive_richardson(n)

function [s A B] = internal(n)
	% set up the system
	 A = (sym(2).^-(0:n-1)'*ones(1,n)).^(ones(n,1)*2*(0:n-1));
	 B = [eye(n) (sym(2).^-(0:n-1)').^(2*n)];
	 % solve
	 s = A \ B;
end

% tabular
syms S R S_num
if (0 == nargin)
 for n=2:6
  [s A B] = internal(n);
  [num den] = numden(s(1,1:end-1));
  D(n-1,1) = den(1); % = max
  S(n-1,1:size(s,2)-1) = s(1,1:end-1);
%  S_num(n-1,1:size(s,2)-1) = D(n-1,1)*num;
  S_num(n-1,1:size(s,2)-1) = D(n-1,1)*num./den;
  R(n-1,1) = s(1,end);
 end
 A
 B
 [D S_num R]
 [S R]
num
den
else
 s = internal(n);
 coeff = double(s(1,1:end-1));
 residual = double(s(1,end));
end

%{
syms h c fh fch fdh fs c d e r

A= ...
[1 -h^2    ;
 1 -c^2*h^2]
B =  [1 0 h^4;
      0 1 c^4*h^4]
F = [fh; fch; r];

s = A \ B; %(b*F);
s = s(1,:);
subs(subs(s,c,1/2),d,1/4)
pause

syms h c fh fch fdh fs c d e

A= ...
[1 -h^2     -h^4;
 1 -c^2*h^2 -c^4*h^4;
 1 -d^2*h^2 -d^4*h^4]
B =  [1 0 0 h^6;
      0 1 0 c^6*h^6;
      0 0 1 d^6*h^6]
F = [fh; fch; fdh; e];

s = A \ B; %(b*F);
s = s(1,:);
subs(subs(s,c,1/2),d,1/4)
pause
%(e*h^6)/64 - (20*fch)/45 + (64*fdh)/45 + fh/45

syms h fh fch fdh feh fs c d e r
A= ...
[1 -h^2     -h^4 -h^6;
 1 -c^2*h^2 -c^4*h^4 -c^6*h^6;
 1 -d^2*h^2 -d^4*h^4 -d^6*h^6;
 1 -e^2*h^2 -e^4*h^4 -e^6*h^6]
B =  [1 0 0 0 h^8;
      0 1 0 0 c^8*h^8;
      0 0 1 0 d^8*h^8;
      0 0 0 1 e^8*h^8]
%F = [fh; fch; fdh; feh; r]

s = A \ B; %(B*F);
s = s(1,:);
subs(subs(subs(s,c,1/2),d,1/4),e,1/8)
pause

%- (r*h^8)/4096 + (84*fch - 1344*fdh + 4096*feh - fh)/2835

%}
end % function derive_richardson

