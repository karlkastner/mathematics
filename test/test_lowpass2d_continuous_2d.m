
L = 10;
n = L^2;

% for idx=1:numel(fr)
%	Sr(idx) = integral(@(r) besselj(0,(2*pi)*r*fr(idx)).*r.*exp(-a*r),0,inf);
% end


%% exp(-a r) = 

syms f a

%symfun = str2sym('2*pi*a*(a.^2 + 4*pi.^2.*f.^2).^(-3/2)');
%fun(f,a)
%I = int(symfun,f,0,inf)

a = 2;
%fun = @(f,a) 2*pi*a*(a.^2 + 4*pi.^2.*f.^2).^(-3/2);
fun = @(f,a) ((2*pi)*a^2)*(a.^2 + 4*pi.^2.*f.^2).^(-3/2);

fr = fourier_axis(L,n);
fr = fr(fr>=0);
fr = mid(fr);
S = lowpass2d_continuous_pdf(fr,a,1,0);
S(:,2) = fun(fr,a);

plot(fr,S)
quad(@(x) fun(x,a),0,L)
sum(S)*(fr(2)-fr(1))
