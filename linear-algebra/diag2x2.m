% So 7. Feb 00:03:10 CET 2016
%
%% diagonal of stacked 2x2 matrices
function A = diag2x2(D)
	A = zeros(2,2,size(D,2));
	A(1,1,:) = D(1,:);
	A(2,2,:) = D(2,:);
end
