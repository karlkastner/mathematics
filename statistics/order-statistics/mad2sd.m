% 2016-06-03 16:01:13.517443048 +0200
% Karl Kastner, Berlin
%
%% transform median absolute deviation to standard deviation
%% for normal distributed values
%
function sd = mad2sd(mad)
	sd = mad/norminv(3/4);
end

