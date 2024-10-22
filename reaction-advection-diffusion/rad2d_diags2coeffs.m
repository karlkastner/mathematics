% Tue 17 Sep 13:26:28 CEST 2024
% Karl Kastner, Berlin
%
%%  recover coefficients of rad-system from diagonals of a 5-point kernel
% ex dz/dx^2 + ey dz/y^2 + ax dz/dx^2
% [0,1,0,-1,1]	  [dx]	[up
% [1,0,-1,0,1]	  [dy]	 left
% [2,2,0,0,1 ]	= [ax]	 centre
% [0,1,0,1,1 ]	  [ay]	 right
% [1,0,1,0,1 ]	  [c]	 down]
function [c,A] = rad2d_diags2coefficients(diagonals)
s=-1;
A = [ 0, 1, 0,-s/2, 0;
      1, 0,-s/2, 0, 0;
     -2,-2, 0, 0, 1;
      1, 0, s/2, 0, 0;
      0, 1, 0, s/2, 0];
c = (diagonals / A');

end

