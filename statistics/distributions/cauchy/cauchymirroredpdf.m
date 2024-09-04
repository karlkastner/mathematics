% 2024-05-24 17:28:18.469718584 +0200
function p = cauchymirroredpdf(f,f0,s)
	p = (s.*(f.^2 + f0.^2 + s.^2))./(pi*((f.^2 + s.^2).^2 + f0^2.*(f0.^2 - 2*f.^2 + 2*s.^2)));
end

