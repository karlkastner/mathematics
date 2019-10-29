% Sun Feb 26 19:08:56 MSK 2012
% Karl Kästner, Berlin

function g = grad_3d(T, P, u)
	g = zeros(size(T,1),3);
	for idx=1:size(T,1)
		M = [1 P(T(idx,1),:);
		     1 P(T(idx,2),:);
		     1 P(T(idx,3),:);
		     1 P(T(idx,4),:)];
		% calculate basis functions
		iM = inv(M);
		g(idx)=	iM(:,2:end)' *[ u(T(idx,1));
			                u(T(idx,2));
			                u(T(idx,3));
			                u(T(idx,4)) ];
	end % for idx
end % grad_3d

