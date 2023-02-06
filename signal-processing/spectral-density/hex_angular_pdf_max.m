function fmax = hex_angular_pdf(k)
	fmax = 0;
	for idx=0:5
		fmax = fmax+misespdf(2*pi*idx/6,0,k);
	end
	fmax = fmax/6;
end

