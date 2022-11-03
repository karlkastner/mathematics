% 2015-03-30 12:32:22.127024358 +0200
% Karl Kastner, Berlin

nf = 100;
X = (1:100)';
%Y = ones(100,1);
Y = (1:100)';
Yf = trifilt1(Y,nf);
Yf(:,2) = trifilt1(Y,nf+1);
plot(X,[Y, Yf, meanfilt1(Y,nf)]);

if (0)

figure(1);
clf
m = 400;
L=16;
nf= m*(2.^(-4:1/16:4));
N = zeros(size(nf),2);
	for idx=1:length(nf);
		x = linspace(0,L,m*L)';
		y = sin(2*pi*x);
		ym = meanfilt1(y,nf(idx));
		yt = trifilt1(y,nf(idx));
		if (mod(idx,8) == 0)
			subplot(3,7,idx/8);
			plot(x,[ym yt]);
			title(nf(idx)/m);
			ylim([-1 1]);
			xlim([7 9]);
		end % if
		N(idx,:) = [norm(ym(700:900)) norm(yt(700:900))];
	end % for
figure(2);
clf
loglog(nf/m,N)

end

