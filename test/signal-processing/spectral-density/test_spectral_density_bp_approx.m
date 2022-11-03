 r = 1e-3; L=10; n = 1e3; nn=1e3; f0=2; f = linspace(0,10*f0,nn)'; S=spectral_density_bp_approx(f,f0,L,n,1,'r',true); S(:,2) = spectral_density_bp(f,f0,L,n,1,'f',true); clf; plot(f,S(:,1)); hold on; plot(f,S(:,2),'--');

