% Fri 10 Dec 11:22:33 CET 2021
% Karl Kästner, Berlin
%
%% filter with a triangular window along both dimensions
% function z = trifilt2(z,nf)
function z = trifilt2(z,nf)
	if (length(nf)==1)
		nf(2)=nf(1);
	end

	z = padd2(z,nf);
	%z=medfilt2(z,[3,3]);
	z = trifilt1(z,nf(1));
	z = trifilt1(z',nf(2))';
	z=z(nf(1)+1:end-nf(1),nf(2)+1:end-nf(2));
end
