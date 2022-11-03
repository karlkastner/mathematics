% 2022-06-13 10:18:31.647787389 +0200

m = 4*1e3;
n = 25*[1,1];
L = 1*[1,1]
%L = [1,0.25];


%[y,x1,x2] = brownian_motion_2d_laplacian([n,1],L);
%y = brownian_motion_2d_fourier([n,m]);
bridge = false;
y = brownian_motion_2d_fft([n,m],L,bridge);
	
x1 = (0:n(1)-1)'/n(1);
x2 = (0:n(2)-1)'/n(2);

%y = y - y((n(1)/2-1)*n(2) + n(2)/2,:);
%y = y - y(1,1,:);
k1 = 1; %n(1)/2;
k2 = 1; %n(2)/2;
% subtract corner value
y = y - y(k1,k2,:);
s = std(y,[],3);

y_ = y; y_ = y_-y_(1,end,:);
s_ = std(y_,[],3);
s_ = rot90(s_);

% end piece
y_ = y; y_ = y_-y_(end,1,:);
se = std(y_,[],3);


%s(n(1)/2)
%y = reshape(y,n);

%x_ = (1:n(1))'/n(1);
if (0)
%	s_ = L(1)*sqrt(x.*(1-x));
else
%	s_ = sqrt(x1);
end

clf
subplot(2,3,1)
plot(x1,s(:,k2)); %/s(n(1)/2,1))
hold on
plot(x1,se(end:-1:1,1));
%plot(x1,s_(:,k2)); %/s(n(1)/2,1))
%set(gca,'colororderindex',1)
plot(x1,sqrt(x1),'--');

subplot(2,3,2)
%set(gca,'colororderindex',2)
plot(x2,s(k1,:)); %/s(1,n(2)/2))
hold on
plot(x1,se(end,:));
%plot(x2,s_(k1,:)); %/s(1,n(2)/2))
plot(x1,sqrt(x2),'--');
%plot(x,diag(s)); %/s(n(1)/2,n(2)/2))
%plot(x1,s_)
%legend('sd_x','sd_y','sd_r');

if (n(1)==n(2))
subplot(2,3,3)
r = hypot(x1,x2);
plot(r,diag(s))
hold on
plot(r,diag(s_))
plot(r,sqrt(r))
legend('d1','d2','expected')
%imagesc(s)
%axis square
%title('Std');
%colorbar
end

subplot(2,3,4)
r = hypot(x1,x2');
%imagesc(s)
plot(r(:),s(:),'.')
%axis square
%title('Std');
colorbar

subplot(2,3,5)
imagesc(mean(y,3))
axis square
title('Mean');
colorbar

if (0)
subplot(2,3,6)
imagesc(y(:,:,1));
title('BM one instance');
colorbar
end

if (0)

% note : interestingly, repredicting the pattern from the spectral density creates a more naturally looking pattern
subplot(2,4,6)
f=100;
n=1e3;
% TODO, with aspect ratio not 1:1, this does not seem to work properly
L = [1,0.5];
 y = brownian_motion_2d_fft([1e3,1e3]);
  x=linspace(0,1,n)';
 b=cos(2*pi*f*(x+0.025*y));
 imagesc(b>0);
 axis square;
 colormap gray;
% daspect([1,4,1]);
 axis([0,1,0,1]*n(1))
S = periodogram_2d(b-mean(b(:)),L);
subplot(2,4,7);
S = trifilt2(fftshift(S),21);
imagesc(S)
S  = fftshift(S);
subplot(2,4,8);
cla
plot(S(:,1)/sum(S(:,1)))
hold on
plot(S(1,:)/sum(S(1,:)))

end

