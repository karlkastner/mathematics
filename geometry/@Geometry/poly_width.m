% Thu 26 May 09:33:47 CEST 2016

%% width of polygon width holes by surface normals
%% holes / islands separated with NaN
%% order of points of outer boundary must be cw
%% order of points of holes must be ccw
%% note that this function does not give the true width for expanding sections
%% use voronoi polygons for this
function W = poly_width(X,Y)
	X = cvec(X);
	Y = cvec(Y);
	np = length(X);

	% TODO split for islands

	% get direction at polygon points by finite differences
	%dXY = [cdiff(X,'circ'), cdiff(Y,'circ')];
	dXY = [cdiff(X), cdiff(Y)];
	% normalisation not necessary
	dXY = bsxfun(@times,dXY,1./hypot(dXY(:,1),dXY(:,2)));
	
	% get orthogonal (normal) direction
	oXY = [dXY(:,2),-dXY(:,1)];

	% unit vectors pointing inward from point p
	uX = [X, X+oXY(:,1)]';
	uY = [Y, Y+oXY(:,2)]';

%figure()
%	plot(uX,uY,'.-')
%pause

	% get edge coordinates
	eX = [X, [X(2:end); X(1)]]';
	eY = [Y, [Y(2:end); Y(1)]]';
	
	% width at every point
	W = NaN(np,1);

%figure(3)
%clf
	% for each point
	timer = tic();
	tlast = 0;
	parfor idx=1:np
%		if (tic(timer)>tlast+10)
%			tlast = tic();
%			printf('%f %f\n',tlast,idx/np);
%		end
%		idx/np
		% get intersection of lines through unit vector and all segments
		%[void s t] = Geometry.lineintersect(uX(:,idx),uY(:,idx),eX,eY);
		[void s t void2 q den] = Geometry.lineintersect( ...
					 [uX(1,idx); uY(1,idx)], ...
					 [uX(2,idx); uY(2,idx)], ...
					 [eX(1,:);   eY(1,:)], ...
					 [eX(2,:);   eY(2,:)]);

		% check that it cuts the second segment (convexity)
		fdx = (t>0) & (t<1);
%		fdx = (s>0) & (s<1);
%		fdx = (s<0) & (s>1);

		% check that it is pointing away from the boundary
%		fdx = fdx & (abs(s) > sqrt(eps));
%		s   = s.*den;
		fdx = fdx & s > 0;
		fdx = find(fdx);
		W_ = hypot(q(1,:)-uX(1,idx), q(2,:)-uY(1,idx));
		[W_ mdx] = min(abs(W_(fdx)));
%		plot([q(1,fdx(mdx)) uX(1,idx)], [q(2,fdx(mdx)) uY(1,idx)]);
%			text(uX(1,idx),uY(1,idx),num2str(W_));
%		hold on
%		if (0==mod(idx,100))
%			axis equal
%			pause
%		end
%	pause(1)
		% write width
		if (~isempty(W_))
			W(idx) = min(W(idx),W_); %abs(s(fdx)));
		end
%		figure()
		

	end
end % poly_width

