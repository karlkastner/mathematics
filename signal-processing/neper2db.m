% Wed 11 Jul 14:21:36 CEST 2018
%% convert neper to db
function db = neper2db(np)
	% 20/log(10) = log10(exp(2)) ~ 0.87
	db = 20/log(10)*np;
end
