% 2024-05-24 17:28:18.469718584 +0200
function p = cauchywrappedpdf(x,x0,s)
	p = cauchypdf(x,x0,s) + cauchypdf(x,-x0,s);
end

