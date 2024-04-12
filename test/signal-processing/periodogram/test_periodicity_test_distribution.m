% 2022-05-18 17:56:44.292546603
% Shat/Sf ~ nf^d*betarnd(a,b)

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

