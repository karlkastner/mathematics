% one sided approximation in the corner

function [dF dF_num dF_den] = derive_bc_one_sided()

A = -[-1 1 0;
      0 1 -1];
B = [-1 1 -1;
      1 1  1];
A = -[ 1 -1 0;
      1  0 -1];
B = [1 1 1;
     2 4 8];
C = inv([1 0 0;
     0 2 0;
     0 0 6])

A__  = [
    -1     1     0     0     0
     0     1     0     0     0
     0     1    -1     0     0
     0     1     0    -1     0
     0     1     0     0    -1 ];
B__  = [
     -1    1    -1    1
     0     0     0     0
     1     1     1     1
     2     4     8    16
     3     9    27    81
];

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
      1     1      1     1
];
%A = [
%    -1     0     1     0     0
%     0    -1     1     0     0
%     0     0     1     0     0
%     0     0     1    -1     0
%     0     0     1     0    -1 ];
%B = [-2     4    -8    16
%     -1     1    -1     1
%      0     0     0     0
%      1     1     1     1
%      2     4     8    16
%];

C__ = inv( ...
[     1     0     0     0
     0     2     0     0
     0     0     6     0
     0     0     0    24]);

% A F = B C dF
% F = [f(-h) f(0) f(+h) f(2h) f(3h)]
dF = (B*C) \ A
dF_den = 1./min(abs(dF'))'
dF_num = dF .* (dF_den*ones(1,5))

end % function derive_bc_one_sided

