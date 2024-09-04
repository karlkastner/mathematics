% Mon 13 Feb 13:38:01 CET 2023
function S = lognmirroredpdf(f,fm,sf);
	S = 0.5*(lognpdf(f,fm,sf) + lognpdf(-f,fm,sf));
end

