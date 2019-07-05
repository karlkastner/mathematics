% 2016-03-03 18:02:43.008968890 +0100
% Karl Kastner, Berlin
%
%% vandermonde matrix of an integral
%
function row = vanderi_1d(l,r,order)
	l = cvec(l);
	r = cvec(r);
	row = zeros(1,order+1,class(l));
	for idx=1:order+1
		row(:,idx) = 1/idx*(r.^idx - l.^idx);
	end
end

