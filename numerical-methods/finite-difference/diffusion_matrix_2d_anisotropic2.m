% Sun 11 Jul 18:43:53 CEST 2021
% 
function D2 = diffusion_matrix_2d_anisotropic2(n,L,a)
if (0)
	[r,l] = directional_neighbour(a);
	D2      = r + l;
	D2(2,2) = D2(2,2)-2;
end
m=-1;
switch (m)
case {-1}
	D2      = zeros(5);
	D2     = [0;1;2;1;0]*[ 0,1,-2,1,0]
case {0}
	D2      = zeros(5);
	D2(3,:) = [1,0,-2,0,1]
case {1}
	D2      = zeros(5);
	p       = sqrt(2)-1;
	D2(3,:) = [p,1-p,-2,1-p,p]
case {2}
	D2      = zeros(3);
	D2(2,:) = [1,-2,1]
end
	D2 = imrotate(D2,rad2deg(a),'bilinear','crop');
D2	
	D2 = kernel2matrix(n,D2);
end

