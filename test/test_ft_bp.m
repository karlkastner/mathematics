figure(1e3);
 clf;
 fmax=20;
 f=linspace(0,fmax)';
 f0=fmax*0.7;
 plot(f,[ft_bp(f,f0,1,fmax) ft_lp(f,f0,1,fmax) ft_hp(f,f0,1,fmax)],'.');
 vline(bp_max(f0,fmax))
