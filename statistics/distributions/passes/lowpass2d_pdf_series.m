% 2023-03-03 21:02:02.088631189 +0100
% only valid for fr<1
function Sb = lowpass2d_pdf(fr,a,order,m)
	n = m;
	Sb = 0;
	or = 2*pi*fr;
	for l=0:n
		Sb = Sb + (-1)^l*(1/2)^(2*l)/factorial(l).^2*or.^(2*l)*a^-(2*n+2)*factorial(2*l+1);
	end
	Sb = 2*pi*Sb;
end

