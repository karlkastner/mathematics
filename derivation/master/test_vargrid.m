% 2012-03-28 00:22:15
% Karl Kästner, Berlin

N=unique([7 ceil(2.^(3:1/64:12))]);
%N=unique([7 ceil(2.^(3:1:12))]);

Err = [];

for idx=1:length(N)
	n=N(idx);
%	c=1;
%	X = pi*sort(rand(n,1));
%	X = pi*(((1:n)'/(n+1)).^2 + ((1:n)'/(n+1)) + ((1:n)'/(n+1)).^3 + ((1:n)'/(n+1)).^4  + ((1:n)'/(n+1)).^5 + ((1:n)'/(n+1)).^6 + ((1:n)'/(n+1)).^7)/7;
%	X = pi*(1-cos(pi*((1:n)'/(n+1))));
%	X = asin(2*(1:n)'/(n+1)-1) + pi/2; % optimal grid -small initial error
%	pause
	%c=4;
	%X = pi*((exp(c*(1:n)'/(n+1))-1)/(c*exp(1)-1));
	L=2; % pi
%	X = L*(1:n)'/(n+1);
	%X = L*sort(((exp(c*(1:n)'/(n+1))))/exp(c));
%	X = L*sort((exp(c*(1:n)'/(n+1))-1)/(exp(c)-1));
	c=-2;
	X = L*sort(-(exp(-c*(1:n)'/(n+1))-1)/(1-exp(-c)));
%	X = asin(2*(1:n)'/(n+1)-1) + pi/2; % optimal grid -small initial error
%	X = polyval(L*ones(10,1)/10,(1:n)'/(n+1));
%	P = [1/24 1/6 1/2 1 1];
%	P = [1e-12 1/6 1/2 1 1];
%	P = [1e-12 1/6 1/2 1 1];
%	P = [1/6 1/2 1 0];
%	P = [1/2 1 0];
%	P = [1 0];
%	P = [1/720 1/120 1/24 1/6 1/2 1 0];
%	P = [ 1/6 1/2 1 0];
%	P = [  1/2 1 0];
%	P = [  1 0];
%	P = [1/24 1/6 1/2 1 0];
%	X = L*sort(1 - polyval(P,2*(1:n)'/(n+1))/polyval(P,2));
%X
%pause
%	y  = exp(-X);
	y  = exp(-X.^2);
	y2  = (4*X.^2)./exp(X.^2) - 2./exp(X.^2);
	y4  = 12./exp(X.^2) - (48*X.^2)./exp(X.^2) + (16*X.^4)./exp(X.^2);
%	y = 2/(81*sqrt(3))*(27 - 18*X + 2*X.^2).*exp(-X/3);
%	y  = (2-X).*exp(-X/2);
%	y2 = 1./exp(X/2) - (X - 2)./(4*exp(X/2));
%	y4 = 1./(2*exp(X/2)) - (X - 2)./(16*exp(X/2));
%	y= sin(X); y2=-y; y4=y;
%	y=1./X; y2=2./X.^3; y4=24./X.^5;
	D4 = d_vargrid([0; X; pi],4);
%	[L D] = d_vargrid([0; X; pi],2,2);
%	D22 = D*L / D;
	D22 = d_vargrid([0; X; pi],2,2);
	[L ] = d_vargrid([0; X; pi],2,4);
	D24 = L;
	z = D22*y;
	Err(idx,1) = sum(sum((y2(3:end-2) - z(3:end-2)).^2));
	z = D24*y;
	Err(idx,2) = sum(sum((y2(3:end-2) - z(3:end-2)).^2));
	z = D22*(D22*y);
	Err(idx,3) = sum(sum((y4(4:end-3) - z(4:end-3)).^2));
	z = D24*(D24*y);
	Err(idx,4) = sum(sum((y4(5:end-3) - z(5:end-3)).^2));
	z = D4*y;
	Err(idx,5) = sum(sum((y4(4:end-4) - z(4:end-4)).^2));
%	Err(idx,1) = sum(sum((y(3:end-2) + z(3:end-2)).^2));
%	Err(idx,2) = sum(sum((y + D24*y).^2));
%	Err(idx,3) = sum(sum((y - D4*y).^2));
%	Err(idx,3) = sum(sum((y - D2*D2*y).^2));
%	legend('y','D^2_2','D^2_4','D^4')
%	plot(X,y)
if (n==2^8)
	figure(1)
	plot(X,ones(size(X)),'.')
	min(diff(X))
	max(diff(X))
%	pause()
end
if (n==2^10)
	subplot(2,2,1)
	plot(X,[y D22*y D24*y D22*D22*y D24*D24*y D4*y ]);ylim([0 2]); %xlim([0 L])
	subplot(2,2,2)
	plot(X,[y y2 y4]); ylim([0 2]); xlim([0 2])
	subplot(2,2,3)
	plot(X,[y2-D22*y y4-D4*y]); ylim([0 2]); xlim([0 2])
%	pause();
end
end
figure(2)
%Err(:,1)=Err(:,1)./Err(1,1);
%Err(:,2)=Err(:,2)./Err(1,2);
%Err(:,3)=Err(:,3)./Err(1,3);
%Err(:,4)=Err(:,4)./Err(1,4);
%Err(:,5)=Err(:,5)./Err(1,5);
loglog(N,Err);ylim([1e-15 1])
legend('location','southwest','D^2_2','D^2_4','(D_2^2)^2','(D_4^2)^2','D_2^4')
xlabel('number of grid points n')
ylabel('$\Big|\Big|\frac{d^ky_*}{dx^k} \mbox{-} \frac{d^ky_h}{dx^k}\Big|\Big|$','interpreter','latex')
%xlim([1e1 1e3])
ylim([1e-12 1e1])
xlim([3 1e3])
%plot(X)
grid on;
set(gca,'minorgrid','none')
print -depsc rounding_error.eps


