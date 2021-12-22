% 2021-09-24 12:34:13.870661927 +0200
% determine coefficients of the ar2 process from the first two lags of the
% autocorrelation function
function c = ar2_acf2c(r1,r2)
	A = [1,r1;
             r1,1];
	c = A \ -[r1;r2];
end

