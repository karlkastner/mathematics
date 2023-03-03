% Fri  3 Mar 15:01:00 CET 2023
a =0.5;
L=80;
n=20*L;
order = 1;
[x,y,r] = fourier_axis_2d(n/L*[1,1],n*[1,1]);
%r = hypot(x,x');
R=exp(-a*r);
S = ifft2(R);
%imagesc(R);

% scale to 1
S = S/S(1,1);
% bandpass
S = (4*S.*(1-S)).^order;
%S=S(:,1);
%S(:,2)=ifft(R(:,1));
% S=S./sum(S);
% plot(real(S))

%fx = fourier_axis(x);
% Sr = [];
% for idx=1:length(fx);
% Sr(idx,1) = integral(@(x) besselj(0,2*pi*x*fx(idx)).*x.*exp(-abs(x)),0,inf);
%  end;
% S=S./S(1,:);
% plot(fx,S);
fx = fourier_axis(x);
Sr = S(:,1);
Sr(:,2) = bandpass2d_pdf_exact(fx,a,order);
%hold on;
clf
plot(fx,Sr(:,1))
hold on
plot(fx,Sr(:,2),'--');
fc = bandpass2d_pdf_exact_mode(a,order);
vline(fc);
rms(Sr(:,1)-Sr(:,2))./rms(Sr(:,1))
xlim([0,L]/2)
Sc = bandpass2d_pdf_exact_scale(L,n,a,order);
%fr=0:0.01:10; a=10; order=2; S=bandpass2d_pdf_exact(fr,a,order); plot(fr,S);
% note that the old derivation of the 2d bandpass (omitting the bessel function)
% is identical
p = spectral_density_bandpass_continuous_max2par(fc,Sc,1)
S = spectral_density_bandpass_continuous(fx,fc,p)/Sc;
plot(fx,S);
rms(Sr(:,1)-S)./rms(Sr(:,1))

