% 2023-01-29 13:57:27.487917043 +0100
function f = hex_angular_pdf(x,mu,k)
	f = 0;
	for idx=0:5
		f = f+misespdf(x+2*pi*idx/6,mu,k);
	end
	f = f/6;
end

