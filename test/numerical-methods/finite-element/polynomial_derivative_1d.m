% Wed Jul 11 19:46:44 MSK 2012
% Karl Kästner, Berlin

% in : C  = [a0, a1, ..., an]
% out: dC = [a1, 2*a2, ...,n*an, 0]    
%
%
function C_dx = polynomail_derivarive_1d(C)
%	if (nargin()<2)
%		n = 1;
%	end
%	n = size(C,2);
	C_dx = zeros(size(C)); %n,size(C,2));
	for idx=1:size(C,2)-1
		C_dx(:,idx) = (idx)*C(:,idx+1);
	end
end % derivative_1d

