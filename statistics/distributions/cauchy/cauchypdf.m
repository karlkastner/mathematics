% 2024-05-24 17:28:18.469718584 +0200
function p = cauchypdf(x,x0,s)
	p = s./(pi*(s.^2 + (x-x0).^2));
%	p = 1./(pi*(gamma+((x-mu)).^2));
end

