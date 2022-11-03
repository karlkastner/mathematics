% Wed 28 Sep 15:56:53 CEST 2022
nsample = 4e4;
n  = 1e1*[1,1];
nf = 5;
mask = ones(n);
p = 0.05;
if (1)
mask = 0.*mask;
mask(1:10,1:10) = 1;
end
q = approximate_ratio_quantile(mask,nf,nsample,1-p);

imagesc(fftshift(q))
