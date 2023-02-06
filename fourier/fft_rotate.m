% Tue 31 Jan 14:05:43 CET 2023
function S = fft_rotate(S, angle_deg)
	n = size(S);
	n_ = n;
	S = fftshift(S);
	% add one tap for domains with even size, so that mean is centred
	if (mod(n(1),2) == 0)
		S = [S; zeros(1,n(2))];
		n_(1) = n_(1)+1;
	end
	if (mod(n(2),2) == 0)
		S = [S, zeros(n_(1),1)];
	end
        S = imrotate(S,angle_deg,'crop','bilinear');
	% strip the added rows/cols
	if (mod(n(1),2) == 0)
		S = S(1:end-1,:);
	end
	if (mod(n(2),2) == 0)
		S = S(:,1:end-1);
	end
	S = ifftshift(S);
end

