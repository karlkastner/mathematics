% 2022-07-05 10:54:25.047696098 +0200
% Karl Kästner, Berlin
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% b  : binary image of thresholded pattern
% kk : array with number of patches with size identical to index
%      kk[1] : number of patches with size 1
%      kk[2] : number of patches with size 2
%      ...
%
function [pa,pp,np,area,pr,radius] = patch_size_2d(b,L)
	% TODO extend to circular bc
	b = padd2(b,1,0);
	s = size(b);
	visited = false(s);
	n = numel(b);
	kk = zeros(n,1); %0;
	np = 0;
	pp = zeros(n,2);
	stack = zeros(n,2);
	for idx=2:s(1)-1
	 for jdx=2:s(2)-1
		[k,ne] = parse_patch_iterative(idx,jdx);
		if (k~=0)
			if (k>length(kk))
				kk(end+1:k) = 0;
			end
			np = np+1;
			kk(k) = kk(k)+1;
			pp(np,:) = [k,ne];
		end
	 end % for idx
	end % for jdx
	%np = sum(kk);
	pp(np+1:end,:) = [];
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

function [k,ne] = parse_patch_iterative(i,j)
	k = 0;
	ne = 0;
	%e = [];
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
			if (~b(i-1,j)) ne=ne+1; end %e(end+1,1:2) = [i,j,i-1,j]; end 
			stack(tail,:) = [i-1,j];
			tail=tail+1;
			if (~b(i+1,j)) ne=ne+1; end %e(end+1,1:2) = [i,j,i+1,j]; end 
			stack(tail,:) = [i+1,j];
			tail=tail+1;
			if (~b(i,j-1)) ne=ne+1; end % e(end+1,1:2) = [i,j,i,j-1]; end 
			stack(tail,:) = [i,j-1];
			tail=tail+1;
			if (~b(i,j+1)) ne = ne+1; end % e(end+1,1:2) = [i,j,i,j+1]; end 
			stack(tail,:) = [i,j+1];
		end % if
	end % while
end % parse_patch_iterative


function [k,e] = parse_patch_recursive(i,j,d,e)
	try
	if (b(i,j) && ~visited(i,j))
		k = 1;
		visited(i,j) = true;
		% TODO add to edge list
		if (~b(i-1,j)) e(end+1,1:2) = [i,j,i-1,j]; end 
		k = k + parse_patch_recursive(i-1,j,d+1,e);
		if (~b(i+1,j)) e(end+1,1:2) = [i,j,i+1,j]; end 
		k = k + parse_patch_recursive(i+1,j,d+1,e);
		if (~b(i,j-1)) e(end+1,1:2) = [i,j,i,j-1]; end 
		k = k + parse_patch_recursive(  i,j-1,d+1,e);
		if (~b(i,j+1)) e(end+1,1:2) = [i,j,i,j+1]; end 
		k = k + parse_patch_recursive(  i,j+1,d+1,e);
	else
		k = 0;
	end
	catch
		d
		k = 0;
	end
	
end

end

