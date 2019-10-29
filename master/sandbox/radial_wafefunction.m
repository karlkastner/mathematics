2:D
E_n = -1/2(n - 1/2)^2

beta_n = 1/(n-1/2)
R_10 = 2 \beta_1 exp(-\beta_1 \rho)
R_20 = 2 \beta_2 /3^1/2 ( 1 - 2 \beta_2 \rho ) exp(-\beta \rho)
R_30 = \beta_3 / 5^1/2 ( 2  - 8 \beta_3 \rho + 4 \beta_3^2 \rho^2 ) exp(-\beta_3 \rho )

% 2D
beta = @(n) 1/(n-1/2)
f1 = @(r) 2*beta(1)*exp(-beta(1)*r)
f2 = @(r) 2*beta(2)/(sqrt(3))*( 1 - 2*beta(2)*r ).*exp(-beta(2)*r)
f3 = @(r) beta(3) /(sqrt(5))*( 2  - 8*beta(3)*r + 4*beta(3)^2*r.^2 ).*exp(-beta(3)*r )

% 3D
f1 = @(r) 1/sqrt(pi).*exp(-r)
f2 = @(r) 1/(4*sqrt(2*pi))*(2 - r).*exp(-r/2)
f3 = @(r) 1/(8*sqrt(3*pi))*(27-18*r+2*r.^2).*exp(-r/3)

fzero(@(r) f1(r).^2-1e-30, 100)
fzero(@(r) f2(r).^2-1e-30, 100)
fzero(@(r) f3(r).^2-1e-30, 100)

R1 = @(r) r.^2.*f1(r).^2

R=linspace(0,100,100)
plot(R,([f1(R); f2(R); f3(R)])); xlim([0 34]); grid on


