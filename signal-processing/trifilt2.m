% Fri 10 Dec 11:22:33 CET 2021
function z = trifilt2(z,nf)
	z = padd2(z,nf);
	%z=medfilt2(z,[3,3]);
	z = trifilt1(z,nf);
	z = trifilt1(z',nf)';
	z=z(nf+1:end-nf,nf+1:end-nf);
end
