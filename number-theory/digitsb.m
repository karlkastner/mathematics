%  2016-09-22 19:50:28.877439254 +0200
% Karl Kastner, Berlin
%
%% number of digits with respect to specified base
%
function p = digitsb(x,base)
	if (nargin()<2)
		base = 10;
	end
	p = floor(log(abs(x))/log(base))+1;
	p(x==0) = 0; 
	%p = floor(log(abs(x))/log(base))
	%x = x.*(base.^-p);
end

