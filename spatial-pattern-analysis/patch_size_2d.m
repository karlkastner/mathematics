% 2022-07-05 10:54:25.047696098 +0200
%
% b  : binary image of thresholded pattern
% kk : array with number of patches with size identical to index
%      kk[1] : number of patches with size 1
%      kk[2] : number of patches with size 2
%      ...
%
function [pa,area,pr,radius] = patch_size_2d(b,L)
	b = padd2(b,1,0);
	s = size(b);
	visited = false(s);
	n = numel(b);
	kk = zeros(n,1); %0;
	stack = zeros(n,2);
	for idx=2:s(1)-1
	 for jdx=2:s(2)-1
		k = parse_patch_iterative(idx,jdx);
		if (k~=0)
			if (k>length(kk))
				kk(end+1:k) = 0;
			end
			kk(k) = kk(k)+1;
		end
	 end % for idx
	end % for jdx

	pa = kk./sum(kk);
	area   = [1:n];
	radius = 0:ceil(sqrt(n)/pi);
	%  cumulative distribution of area
	Ka = cumsum(kk);	
	Ka = Ka/Ka(end);

	% cumulative distribution of radii
	Kr = interp1(sqrt(area)/pi,Ka,radius,'linear');

	% density and raidus at mid-points
	pr     = diff(Kr);
	radius = mid(radius);

	if (nargin()>1)
		dx = rvec(L)./size(b);
		area = area*dx(1)*dx(2);
		radius = dx*radius;
	end

function k = parse_patch_iterative(i,j)
	k = 0;
	head = 1;
	stack(1,1:2) = [i,j];
	tail = 1;
	while (head <= tail)
		% pop the stack
		[ij] = stack(head,:);
		i = ij(1);
		j = ij(2);
		head = head+1;
		if (b(i,j) && ~visited(i,j))
			k = k+1;
			visited(i,j) = true;
			% push the neighbours
			tail=tail+1;
			stack(tail,:) = [i-1,j];
			tail=tail+1;
			stack(tail,:) = [i+1,j];
			tail=tail+1;
			stack(tail,:) = [i,j-1];
			tail=tail+1;
			stack(tail,:) = [i,j+1];
		end % if
	end % while
end % parse_patch_iterative


function k = parse_patch_recursive(i,j,d)
	try
	if (b(i,j) && ~visited(i,j))
		k = 1;
		visited(i,j) = true;
		k = k + parse_patch_recursive(i-1,j,d+1);
		k = k + parse_patch_recursive(i+1,j,d+1);
		k = k + parse_patch_recursive(  i,j-1,d+1);
		k = k + parse_patch_recursive(  i,j+1,d+1);
	else
		k = 0;
	end
	catch
		d
		k = 0;
	end
	
end

end

