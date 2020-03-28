% Wed Feb 27 00:36:42 MSK 2013
% Karl Kästner, Berlin

function test_circmesh()
	opengl neverselect;
	n = 31;
	[P T B X] = circmesh(n,2.0);
	plot(P(:,1), P(:,2),'.')
	axis equal
	hold on
	patch([ P(T(:,1),1)  P(T(:,2),1)  P(T(:,3),1)]', ...
			[ P(T(:,1),2)  P(T(:,2),2)  P(T(:,3),2)]',ones(size(T,1),1)'); %,[0 0 0',varargin{:});
end


