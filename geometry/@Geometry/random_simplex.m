% Thu 11 Jan 14:10:34 CET 2018
%
%% random point inside of a triangle
% see triangle point picking on wolfram
%
% m : number of vertices (number of dimension + 1)
% n : number of random points inside the simplex
%
function [c P0 R] = random_simplex(n,m,varargin)
%	P1 = varargin{1};
%	n  = length(varargin{1});
%	m  = nargin;
	
	% random coordinates for points 2-n
	R = rand(n,m-1);
	
	% barycentric coefficients
	c = [R,1-sum(R,2)];
	% with transformation: c = [R,1];
	% points where c(3) < can be mirrored at (0.5,0.5,0)
	% TODO the mirroring does not work in 3d
	% -> mirror at two axes/points?
	fdx = c(:,end) < 0;
	nf  = sum(fdx);
	%d          = [3/(m-1)*ones(nf,m-1),zeros(nf,1)]-c(fdx,:);
	c(fdx,:) = [2/(m-1)*ones(nf,m-1),zeros(nf,1)]+c(fdx,:);
%	c(fdx,:) = [2/(m-1)*ones(nf,m-1),zeros(nf,1)]-c(fdx,:);

	% transform to bacycentric coordinates
	if (length(varargin)>0)
		P0 = zeros(n,m-1);
		for idx=1:m
			Pi = varargin{idx};
			P0 = P0 + bsxfun(@times,c(:,idx),Pi);
			%P0 = P0 + bsxfun(@times,R(:,idx),(Pi-P1));
		end
	end
end

