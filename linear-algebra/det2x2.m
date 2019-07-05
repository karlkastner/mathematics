% Mo 21. Dez 00:21:43 CET 2015
% Karl Kastner, Berlin
%% 2x2 matrix inverse of 2x2 matrices stacked along dim 3
function det = det2x2(A)
	det  = (A(1,1,:).*A(2,2,:) - A(1,2,:).*A(2,1,:));
end

