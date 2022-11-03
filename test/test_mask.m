% 2022-09-23 15:54:22.164843714 +0200

% repeat
% generate gaussian noise, mask, compute S, compute Sbar, compute ratio S/Sbar
% compute moments of the ratio


if (0)

n2 = 1e3;

L=10;
n = 1e4;
x = innerspace(0,100,n)';
l = 1;
%0.97);
y = cos(2*pi*x/l);
%y = randn(n,1);
y = y./rms(y);
y(:,2) = y.*(x<L);
% +0*randn(length(x),1).*(x>L);
%y(:,3) = 0;
m = sum(x<L);
y3 = randn(n,n2);
y3 = y3./rms(y3);
y3(1:m,:) = repmat(y(1:m,1),1,n2);
y3 = y3./rms(y3);
S3 = mean(abs(fft(y3,[],1)).^2,2);

%y(1:10:end,3) = y(1:10:end,1);



if (0)
% sellect random subset
id = randperm(n,m);
y(id,3)=y(id,1);
y(:,4)=0;
end
%y(id,1);

y=y./rms(y);

subplot(2,2,1);
semilogy([abs(fft(y,[],1)).^2,S3]);
legend('unmasked','masked','masked+whitened')

y 

subplot(2,2,2);
end

%y    = zeros(n,1);
%y(1) = 1;
% e^2 = sum ei ej = sum si^2 + sum rho_ij si sj
% - correlation matrix for
% plot(y) 

% cij = corr yi yj, where y = |(F w e)|^2
% cij = E[yi yj] = E[ |(F w e)|^2_i |(F w e)|^2_j)
% -> C = |(F w)|^2 * |(F w)^*|^2

% to test correlation matrix
%	define a mask
%	define a number of random input vectors
%	compute fft of masked vectors
%	compute cov of vectors
%	check that this matches the cov computed via ft
% to test dof : std of masked S-values

% -> then check periodicity test with mask

% nt=10; m=1e2; n = 1e4; for idx=1:n; x= randn(m); x(end/2:end,:) = 0; x(:,end/2:end)= 0; p(idx,1) = periodogram_test_periodicity_2d(x,[],nt); end; sum(p<0.05)/n

