% Thu May  3 02:53:01 MSK 2012
% Karl Kästner, Berlin

function [T n] = restore_cw(P, T, B)
	% restore clockwise ordering
	n = 0;
	for idx=1:size(T,1)
		A = [ 	1 P(T(idx,1),:);
			1 P(T(idx,2),:);
			1 P(T(idx,3),:)];
		if (det(A) < 0)
			help = T(idx,1);
			T(idx,1) = T(idx,2);
			T(idx,2) = help;
			n = n+1;
		end
	end
end % restore_cw

