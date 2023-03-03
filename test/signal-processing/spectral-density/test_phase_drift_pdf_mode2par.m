% Tue 18 Jan 10:33:57 CET 2022
function fail = test_spectral_density_brownian_phase_mode2par()
	fc = [1,2,3]';
	
	Sc = [0.5,0.7,1];
	
	f0 = [];
	s  = [];
	for idx=1:length(fc)
	 for jdx=1:length(Sc)
		[f0(idx,jdx), s(idx,jdx)] = spectral_density_brownian_phase_mode2par(fc(idx),Sc(jdx));
	 end
	end
	f0
	s
	[fc_,Sc_] = spectral_density_brownian_phase_mode(f0,s)
	
	res = [(fc_ - fc).^2 + (Sc - Sc_).^2];
	
	fail = max(abs(res(:))) > 1e-4;
end

