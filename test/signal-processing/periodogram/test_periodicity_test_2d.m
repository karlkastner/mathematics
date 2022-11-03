% Shat/Sf ~ nf^d*betarnd(a,b)

if (0)
n=100;
nf = 3;
nf2 = 5;
b = randn(n,n);
Shat = abs(fft(b)).^2;
Shat = fftshift(Shat);

S = Shat;
if (nf > 1)
	S = meanfilt1(S,nf);
end
if (nf2 > 1)
	S = meanfilt1(S',nf2)';
end
S = meanfilt2(Shat,nf);
%if (m==1)
%S = meanfilt1(Shat,nf);
%else
%S = meanfilt2(Shat,nf);
%end

r = Shat./S;
%r = r(1:end/2,:);

nn = numel(r);
p = (1:nn)'/(nn+1);

if (1 == m)
        a = 1;                                                          
        b = (nf-1);       
	q = nf*betainv(p,a,b);
else
        a = 1;                                                          
        b = (nf*nf2-1);       
	q = nf*nf2*betainv(p,a,b);
%         p1 = 1-betacdf(ratio/nf^2,a,b);   
end
	plot(p,[q,sort(flat(r))]);
end

if (1)

nt=1e4;
n = 100;
p =[];
L = 1;
nf = 5;
for idx=1:nt
	b = randn(n,n);
	[p(idx),ratio,maxShat,mdx,fdx,S] = periodogram_test_periodicity_2d( b, L, nf);
end
sum(p>0.05)/nt
sum(p>0.1)/nt
sum(p>0.2)/nt

end

