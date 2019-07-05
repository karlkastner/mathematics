% Sa 28. Nov 12:07:01 CET 2015
% Karl Kastner, Berlin
%
%% edge length
function H = edge_length(E,P)
	H = hypot(P(E(:,1),1) - P(E(:,2),1), ...
		  P(E(:,1),2) - P(E(:,2),2));
end

