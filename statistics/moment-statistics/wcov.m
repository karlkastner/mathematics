% Thu 10 Aug 18:27:42 CEST 2017
%
%% covariance of two vectors when samples are weighted
function v_12 = wcov(w1,x1,w2,x2)
	v_12  =   wmean(w1.*w2,x1.*x2) ...
               -  wmean(w1,x1).*wmean(w2,x2);
end

