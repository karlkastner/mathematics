% Mo 30. Mär 11:52:39 CEST 2015
% Karl Kastner, Berlin
%
%% cumulative sum, setting nan values to zero
%
function x = nancumsum(x)
	fdx = isnan(x);
	x(fdx) = 0;
	x = cumsum(x);	
end

