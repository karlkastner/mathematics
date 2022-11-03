% 2021-05-17 19:20:51.352432728 +0200

% note : padding the same values is better than padding zeros
%	 padding 

namedfigure(1,'oversampling');
clf
e = 1e-15;
e = 0;
f = @(x) sign(x)+e*randn(size(x))
f = @(x) (x)+e*randn(size(x))
clf;
n=16*16;
id  = 1:n/2;
m   = 2;
id_ = 1:m*n/2;
subplot(2,2,1);
L=16;
x=innerspace(0,L,n)';
d=[-1/2,0,1/2];
d = [0];
for idx=1:length(d);
	y=cos(2*pi*(x/L*(L-d(idx))));
	ff = fft(f(y),m*n);
	subplot(2,2,1);
	stairs(abs(ff(id_)),'.','linewidth',3);
        [mv,mdx] = max(abs(ff(id_)));
	set(gca,'xscale','log','yscale','log')     
	hold on;
	subplot(2,2,2);
	stairs(angle(ff(id_)),'.','linewidth',3);
	vline(mdx)
	set(gca,'xscale','log')     
	hold on
end

m=16;
d_ = 0.25;
for idx=1:3;
	y=cos(2*pi*(x/L*(L-d_)));
	switch(idx)
	case {1}
		m_ = 1;
	case {2}
		m_ = m;
	case {3}
		y = repmat(y,m);
		m_ =m;
	end
	ff = fft(f(y),m_*n);
	ff = ff(1:m_*n/2);
	[fa,Ta] = fourier_axis(L*m_,m_*n);
fa(2)
	fa = fa(1:m_*n/2);
	subplot(2,2,3);
	stairs(fa,abs(ff),'.','linewidth',3);
        [mv,mdx] = max(abs(ff));
	set(gca,'xscale','log','yscale','log')     
	hold on;
	subplot(2,2,4);
set(gca,'colororderindex',idx)
	stairs(fa,angle(ff),'.','linewidth',3);
	vline(mdx)
	set(gca,'xscale','log')     
	hold on
end
if (0)
figure(2);
clf
win = tukeywin(n,0.25);
w = win/mean(win);
subplot(2,2,3)
% d=[-1/2,0,1/2];
 for idx=1:2 %length(d);
 y=cos(2*pi*(x/L*(L-0.5)));
if ( 1 == idx)
 ff = ((fft(f(y),m*n)));
else
 ff = ((fft(w.*f(y),m*n)));
end
 subplot(2,2,3);
 stairs(abs(ff(id_)),'.','linewidth',3);
        [mv,mdx] = max(abs(ff(id_)));
 set(gca,'xscale','log','yscale','log')    
 subplot(2,2,4);
 stairs(angle(ff(id_)),'.','linewidth',3);
	vline(mdx)
 hold on;
 set(gca,'xscale','log')     
 end;

% clf; x = innerspace(0,1,1024); x=(1-cos(2*pi*x))/2; for idx=0:3; y=1-(1-x).^idx; loglog(abs(fft(y))); hold on; end
end
