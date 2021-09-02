% Wed 23 Jun 16:32:43 CEST 2021
% note : this seems already bias-corrected
function a = autocorr2(x) %,correct_bias)
	f = fft2(x);
	%f = ifftshift(abs(f_).^2);        
	f = f.*conj(f);                                                     
	a = ifft2(f);
	% bias correction
%	if (correct_bias)
%	n = size(a);
%	if (mod(n(1),2) == 1)
%		s1 = (1:n)
%	else
%		s1 = (n(1)/2)./(n(1)/2-[0:n(1)/2-1,n(1)/2:-1:1]);
%	end
%	if (mod(n(2),2) == 1)
%		s2 = (1:n)
%	else
%		s2 = (n(2))./(n(2)-[0:n(2)/2-1,n(2)/2:-1:1]);
%	end
%	%a = a.*s1'; %.*s2;
%	%a = a.*s2.^1;
%	end

	a = fftshift(a);
	
%        f = ifft(f')';
        %a=f;
        % f=f/max(f(:));
end
