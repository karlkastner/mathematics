% this is just a fist order lowpass filter
function S = ornstein_uhlenbeck_spectral_density(f,theta)
	fc = theta/(2*pi);
	scale = 2/(pi*fc);
	S = scale./(1 + (f./fc).^2);
end

