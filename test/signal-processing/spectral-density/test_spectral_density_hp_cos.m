
 
 fc_pattern = 1;
 fmu = 3*fc_pattern;
 fc_lp = fmu/5; 

 x =linspace(0,10*fc_pattern)';
 fx=fourier_axis(x);
 S = sd_lp_cos(fx,fc_lp,1);
 S2 = sd_lp_cos(fx,fc_lp,2);
% S2= sd_lp_cos(fx,1/sqrt(2));
 plot(fx,[S,S2],'.-');
 vline(fc_pattern)

% plot(x,ifft(S.^0.5),'.-');
% xlim([0,3])

L=100; e=randn(1e3,1); [y,S] = highpass1d_fft_cos(e,1,L,2); plot(S); xlim([0,3]); plot(abs(fft([y,e])).^2);
