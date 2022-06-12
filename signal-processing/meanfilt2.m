% Fri 10 Dec 11:22:33 CET 2021
% Karl Kästner, Berlin
%
%% filter with a rectangular window along both dimensions
% function z = meanfilt2(z,nf)
function z = meanfilt2(z,nf)
	z = padd2(z,nf);
	z = meanfilt1(z,nf);
	z = meanfilt1(z.',nf).';
	z = z(nf+1:end-nf,nf+1:end-nf);
end
