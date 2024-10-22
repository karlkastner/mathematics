% 2023-07-20 10:15:04.929342964 +0200
% Tue  7 Nov 10:22:34 CET 2023
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
% downsampling aka deflation matrix
function A = downsampling_matrix(n,mode)
	switch(mode)
	case{'pairwise'}
		% proper variance scaling, but shifts half a grid point
		i = (1:n)';
		j = flat([1;1]*(1:n/2));
		A = sparse(i,j,0.5,n,n/2);
	case {'multigrid','mg'}
		% improper variance scaling, but keeps every other point without shifting
		w = [1,2,1]/4;
		% w = [1/6,2/3,1/6];
		A = spdiags(ones(n,1)*w,-1:1,n,n);
		A = A(1:2:end,:);
		A(1,end) = 0.25;
		%buf = [id,2*id-1,ones(n/2,1); id, 2*id, 2*ones(n/2,1); id, mod(2*id+1,n),ones(n/2,1)];
		%P = sparse(buf(:,1),buf(:,2),buf(:,3));
	otherwise 
		error('here');
	end
end

