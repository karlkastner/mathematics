% Sun Jul 15 18:36:14 MSK 2012
% Karl Kästner, Berlin
%
%% fast fourier transform for complex input data
%%
%% input:
%% F : data in real space
%%
%% output :
%%
%% F : fourier transformation of F
%%
function F = fft_man(F)
	n = size(F,1);
	if(n < 2)
		F = F;
		return
	end
	if (bitand(n,n-1) & n)
		error('fft_man','Input vector lenght must be a power of two');
	end
	m = n/2;

	F_even  = fft_man(F(2:2:n,1));
	F_odd   = fft_man(F(1:2:n-1,1));
	% in case of NaN, one can assume even == odd, but
	% then samples are differently weighed
	% if all but one samples are nan in one half, this sample will give weight to the mean
	% as much as all other samples together
	% alternatively, put to zero
%	flag = isnan(F_even);
%	F_even(flag) = -conj(F_odd(flag));
%	F_even(flag) = conj(F_odd(flag));
%	F_even(flag) = 0;
%	flag = isnan(F_odd);
%	F_odd(flag) = -conj(F_even(flag));
%	F_odd(flag) = conj(F_even(flag));
%	F_odd(flag) = 0;

	% tricky the transpose here, use
	% either -i together with .'
	%     or +i together with '
	D = exp(-1i*2*pi/n).^(0:m-1).';
	F   = [D.*F_even + F_odd;
              -D.*F_even + F_odd];
	k = m/2;
end % fft_man

