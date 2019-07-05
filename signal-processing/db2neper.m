% Wed 11 Jul 14:21:42 CEST 2018
%% convert decibel to neper
function np = db2neper(db)
	% log(10)/20 = 1/log10(exp(2)) ~ 0.12
	np = log(10)/20*db;
end

