% Mo 21. Dez 00:21:43 CET 2015
% Karl Kastner, Berlin
%% 2x2 inverse of stacked matrices
function Ai = inv2x2(A)
	if (~issym(A))
		Ai = zeros(size(A));
	end
%	a = A(1,1,:);
%	b = A(1,2,:);
%	c = A(2,1,:);
%	d = A(2,2,:);

	idet = 1./det2x2(A);
%(a.*d - b.*c);
%	Ai        = zeros(size(A));
	Ai(1,1,:) =  idet.*A(2,2,:);
	Ai(1,2,:) = -idet.*A(1,2,:);
	Ai(2,1,:) = -idet.*A(2,1,:);
	Ai(2,2,:) =  idet.*A(1,1,:);
end % inv2x2

