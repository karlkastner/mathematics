% 2025-03-11 14:50:15.776691865 +0100
%
% approximate the semi-discrete density S_L1 with finite domain size L1,
% by averaging the infinite-dimensional density S_inf with infinite domain L-> inf
% by averaging the semi-dicrete density in the domain with size m*L1
%
%	S_L1 = 1/(\Delta{f}) \int_{(i-1/2)\Delta{f}}^{(i+1/2)\Delta{f}} S_inf df
%	midpoint rule with 1-based matlab style array indexing:
%           ~ 1/m sum_(j=1)^m S_mL(i*m + j - (m+1)/2)
%	df = 1/L1
%
%	m has to be odd, for m = 1 no oversampling takes place
%
% TODO sample beyond last tap should not be circularly wrapped but considered 0
function S_L = resample_density_by_averaging(S_mL,m,dim)
	if (mod(m,2) ~= 1)
		error('m has to be odd');
	end
	if (nargin()<3)
		dim = 1;
	end
	S_  = 0;
	for idx=1:m
		S_  = S_ + circshift(S_mL,idx-(m+1)/2,dim);
	end
	if (1 == dim)
		S_ = S_(1:m:end,:);
	else
		S_ = S_(:,1:m:end);
	end
	S_L = S_/m;
end

