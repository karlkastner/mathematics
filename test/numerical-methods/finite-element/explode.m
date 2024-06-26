% Sun Jul 29 17:51:19 MSK 2012
% Karl Kästner, Berlin

function [P_ T_] = explode(P, T, x0, c)
	if (isempty(x0))
		x0 = mean(P);
	end
	%P_ = P(T(idx,:))
	P_ = [];
	T_ = [];
	for idx=1:length(T)
		m = mean(P(T(idx,:),:))
		P(T(idx,:),:)
	(c-1)*ones(4,1)*m 
		P_ = [P_; P(T(idx,:),:) + (c-1)*ones(4,1)*(m-x0) ];
		T_ = [T_; (idx-1)*4 + [1 2 3 4]];
	end
%	x0_ = ones(size(P,1),1)*x0;
%	P = c*(P - x0_) + x0_;
end

