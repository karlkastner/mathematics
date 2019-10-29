% Tue May  1 01:55:50 MSK 2012
% Karl Kästner, Berlin

function display_2d(mesh, t, c, pdx, varargin)
	opengl neverselect;
	P = mesh.P;
	T = mesh.T;
	B = mesh.Bc;

	h = ishold();
	hold on;
	if (nargin < 6 || isempty(pdx))
		pdx=(1:size(T,1))';
	end
	if (nargin < 5 || isempty(c))
		c = 0.5*ones(length(pdx),1);
	end

	% triangles
	patch([ P(T(pdx,1),1)  P(T(pdx,2),1)  P(T(pdx,3),1)]', ...
			[ P(T(pdx,1),2)  P(T(pdx,2),2)  P(T(pdx,3),2)]',c',varargin{:});

	% boundaries
	if (nargin > 3 && bitand(t,8) == 8)	
		plot([P(B(:,1),1) P(B(:,2),1)]', [P(B(:,1),2) P(B(:,2),2)]','o-r','Linewidth',2);
	end
	% point label
	if (nargin > 3 && bitand(t,1) == 1)	
		for idx=1:size(P,1)
			text(P(idx,1), P(idx,2), num2str(idx))
		end
	end
	% triangle label
	if (nargin > 3 && bitand(t,2) == 2)	
		for idx=1:length(pdx)
			text(mean(P(T(pdx(idx),:),1)),mean(P(T(pdx(idx),:),2)), num2str(pdx(idx)),'color',[0 0 1])
		end
	end
	% boundary label
	if (nargin > 3 && bitand(t,4) == 4)	
		for idx=1:size(B,1)
			text(mean(P(B(idx,1:2),1)),mean(P(B(idx,1:2),2)), num2str(idx),'color',[1 0 0])
		end
		set(gca,'xtick',[])
		set(gca,'ytick',[])
	end

	if (~h)
		hold off;
	end
end % display_2d

