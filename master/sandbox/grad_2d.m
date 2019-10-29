% Sun Feb 26 19:08:36 MSK 2012
% Karl Kästner, Berlin

function g = grad_2d(T, P, u)
	g = zeros(size(T,1),2);
	for idx=1:size(T,1)
		M = [1 P(T(idx,1),:);
		     1 P(T(idx,2),:);
		     1 P(T(idx,3),:)];
		% calculate basis functions
		iM = inv(M);
		g(idx)=	iM(:,2:end)' *[ u(T(idx,1));
			                u(T(idx,2));
			                u(T(idx,3)) ];
	end
end % grad_2d

