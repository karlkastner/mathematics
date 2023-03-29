% 2023-03-13 16:22:14.907808880 +0100
L = 10;
n = 2000;
if (1)
n = 128;
m = 1000;
s21 = 0;
mu21 =0;
s22 = 0;
mu22 =0;
for idx=1:m	
	idx/m
	rng(idx+2*m)
	b_ = brownian_field(0.5,n);
	res1_ = b_ - b_(:,1);
	s21   = s21 + res1_.^2/m;
	mu21  = mu21 + res1_/m;
	res1_ = b_ - b_(1,:);
	s22   = s22 + res1_.^2/m;
	mu22  = mu22 + res1_/m;
end
end
sc = sqrt(1);
b = sqrt(L)*b_;
clf
subplot(2,2,2)
s2 = var(b-b(:,1),[],1)';
s2(:,2) =rms(b-b(:,1),1).^2;
plot(s2);

subplot(2,2,1)
s1 = var(b-b(1,:),[],2);
s1(:,2) =rms(b-b(1,:),2).^2;
plot(s1);

subplot(2,2,3);
plot(mean(sc*s21,1));
hold on
plot(mean(mu21,1));

subplot(2,2,4);
plot(mean(sc*s22,2));
hold on
plot(mean(mu22,2));

