% Mon 13 Feb 13:38:01 CET 2023
function S = normpdf_wrapped(f,fm,sf);
	S = normpdf(f,fm,sf) + normpdf(f,-fm,sf);
end

