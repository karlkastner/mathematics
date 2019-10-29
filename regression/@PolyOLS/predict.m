% Fri Feb 20 13:05:45 CET 2015
% Karl Kastner, Berlin
%% predict polynomial function values
% TODO, weighing
function [Y, Sc, Sp, obj] = predict(obj,X)
	X = (X - obj.x0)/obj.s;
	[Y, A] = obj.predict_(obj.param,X);
%	S2 = diag((A*obj.C)*A');
	% most listerature only treats the trivial case of slope and intercept (SLS)
	% Prediciton error of OLS can be found in bingman
	% (Note his confusing notation to rotate the prediction matrix)
	nsigma = 1;
	c  = tinv(normcdf(nsigma),obj.nsample-obj.nparam);
	S2 = sum((A*obj.C0).*A,2);
	% confindence interval (mean response)
	Sc  = c*sqrt(S2)*obj.serr;
	% predicition interval (individual response)
	Sp = c*sqrt(1+S2)*obj.serr;
%	L  = Y - S;
%	U  = Y + S;
end % predict

