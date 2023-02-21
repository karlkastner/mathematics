% Tue 31 Jan 14:05:43 CET 2023
% Karl Kästner, Berlin
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
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

