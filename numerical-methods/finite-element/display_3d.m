% Thu Jul 26 18:15:15 MSK 2012
% Karl Kästner, Berlin


function display_3d(P,T,Bc,flag,alpha)
	if (nargin() < 5 || isempty(alpha))
		alpha = 0.125;
	end

	hstat = ishold();

	% point label
	if (bitand(flag,8))	
		for idx=1:size(P,1)
			text(P(idx,1), P(idx,2), P(idx,3), num2str(idx),'FontSize',24); hold on
		end
	end


	% plot the tetrahedra
	if (bitand(flag,1))
		if (16 ~= bitand(flag,16))
			h = tetramesh(T,P);
			set(h,'facealpha', alpha);
			hold on;
		else
		for jdx=1:size(T,1);
			h = tetramesh(T(jdx,:), P, 'r');
			set(h,'facealpha', alpha);
			hold on
			axis([0 1 0 1 0 1])
			pause
		end   
		end
	end

	% plot the boundary triangles
	if (bitand(flag,2))
		alpha = 1;
		C = (1:size(Bc,1));
		h = patch( [P(Bc(:,1), 1) P(Bc(:,2), 1) P(Bc(:,3), 1) ]', ...
		       [P(Bc(:,1), 2) P(Bc(:,2), 2) P(Bc(:,3), 2) ]', ...
		       [P(Bc(:,1), 3) P(Bc(:,2), 3) P(Bc(:,3), 3) ]', C );
		set(h,'facealpha', alpha);
%		for idx=1:size(Bc,1)
%			patch( [P(Bc(idx,1), 1) P(Bc(idx,2), 1) P(Bc(idx,3), 1) ]', ...
%			       [P(Bc(idx,1), 2) P(Bc(idx,2), 2) P(Bc(idx,3), 2) ]', ...
%			       [P(Bc(idx,1), 3) P(Bc(idx,2), 3) P(Bc(idx,3), 3) ]', (idx)/size(Bc,1)*[1 1 1] );
%			hold on;	
%			pause
%		end
		hold on;
	end

	% plot points
	if (bitand(flag,4))
		plot3(P(:,1),P(:,2),P(:,3),'*y');
		hold on;
	end

	% hexahedra lable
	if (bitand(flag,16))
		pdx=(1:size(T,1))';
		for idx=1:length(pdx)
			text(mean(P(T(pdx(idx),:),1)), 	mean(P(T(pdx(idx),:),2)), mean(P(T(pdx(idx),:),3)), num2str(pdx(idx)),'color',[1 0 0],'FontSize',24);
		end
	end

%	if (bc_plot)
%		for idx=1:size(Bc,1)
%			fill3(P(Bc(idx,:),1), P(Bc(idx,:),2), P(Bc(idx,:),3), 'g');
%			hold on;
%			axis([0 1 0 1 0 1]);
%		end
%	end

	% plot individual tetras with pause

	if (~hstat)
		hold off;
	end
	grid on;
%	axis([0 1 0 1 0 1])
end % display_3d

