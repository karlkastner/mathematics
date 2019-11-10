% Thu 14 Jun 16:43:54 CEST 2018
function ternary_diagram(pq,label,varargin)
	XY = [-1,0;
               1,0;
               0, 2*sqrt(2/3)];
	trimesh([1,2,3],XY(:,1),XY(:,2),varargin{:});
	hold on
	if (1)
		% plot separator lines
		xyc = barycentric2cartesian(XY,[1, 1, 1]/3);
	%	plot(xyc(1),xyc(2),'*');
		xy0 = barycentric2cartesian(XY,[1, 0, 1]/2);
		plot([xyc(1),xy0(1)],[xyc(2),xy0(2)],'-k',varargin{:});
		plot([xyc(1),XY(2,1)],[xyc(2),XY(2,2)],'--k',varargin{:});
		xy0 = barycentric2cartesian(XY,[0, 1, 1]/2);
		plot([xyc(1),xy0(1)],[xyc(2),xy0(2)],'-k',varargin{:});
		plot([xyc(1),XY(1,1)],[xyc(2),XY(1,2)],'--k',varargin{:});
		xy0 = barycentric2cartesian(XY,[1, 1, 0]/2);
		plot([xyc(1),xy0(1)],[xyc(2),xy0(2)],'-k',varargin{:});
		plot([xyc(1),XY(3,1)],[xyc(2),XY(3,2)],'--k',varargin{:});
	end
	pq = bsxfun(@times,pq,1./sum(pq,2));
	xy0 = barycentric2cartesian(XY,pq).';
	scatter(xy0(:,1),xy0(:,2),'ro','filled');
	axis equal
	axis off
	if (nargin()>1 & ~isempty(label))
		for idx=1:3
			p = 0.5;
			XY_ = bsxfun(@plus,(1-p)*xyc',p*XY);
			text(XY_(idx,1),XY_(idx,2),label{idx});
		end
	end
end

