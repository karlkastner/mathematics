% 2023-03-16 15:11:19.933303369 +0100

n = 16;
m = 1e3;

sx = 1;
sy = 0.5;

s21 = 0;
s22 = 0;
for idx=1:m
	b = brownian_field_scaled(0.5,n,[sx,sy]);
	s22 = s22 + mean((b-b(1,:)).^2,2)/m;
	s21 = s21 + mean((b-b(:,1)).^2,1)'/m;
end

plot(sqrt([s21,s22]));

