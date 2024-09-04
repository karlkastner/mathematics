% Mon 13 Feb 13:38:01 CET 2023
function S = normalmirrored(f,fm,sf);
	S = 0.5*(normpdf(f,fm,sf) + normpdf(f,-fm,sf));
end

