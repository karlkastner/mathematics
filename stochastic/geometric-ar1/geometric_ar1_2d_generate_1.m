% Wed  1 Mar 11:36:59 CET 2023
% Karl Kastner, Berlin
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
%% realization of the spatial geometric ornstein (geometric ar1) process
%% averaged over grid cells
function [a,Ry] = generate_geometric_ornstein_uhlenbeck_averaged(mu_z,sd_z,theta,L,n,m)
	x  = ((0:(n*m)-1)'/(n*m)-0.5)*L;
	x  = fftshift(x);
	% autocorrelation of z = log(a)
	R  = exp(-theta*hypot(cvec(x),rvec(x)));
	% spectral density of z = log(a)
	S  = real(fft2(R));
	% normalize
	df = 1/L;
	S  = S/(sum(S,'all')*df^2);
	% transfer function
	T  = sqrt(S);
%	if (acflag)
%		% this is for testing purposes
%		% note that the acf of exp(z) is not exp(acf(z))
%		e0 = zeros(n*m);
%		e0(1,1) = 1;
%	else
	% white noise
	e0 = randn(n*m);
	% lowpass filter to simulate 2d-ar1 process (ornstein uhlenbeck)
	e0 = real(ifft(T.*fft(e0)));
%	end
	% make log-normal
	a0 = exp(mu + sd_z*e0);
	% allocate output
	a  = zeros(n);
	% average over columns
	for idx=1:m
		% average over rows
		for jdx=1:m
	 		a = a + a0(idx:m:end,idx:m:end);
		 end % for jdx
	end % for idx
	% normalize
	a = a/m^2;

if (0)
% TODO wrap around!
% TODO this does not mimick a radillay symmetric filter,
%      as it is identical to filtering along x and subsequenty along y,
%      leading to a diamond shaped density instead of a circular one
%	-> this can only work when the covariance is gaussian

	% allocate output
	a  = zeros(n);
	x  = ((0:n-1)'/n-0.5)*L;
	x  = fftshift(x);
	y  = x;


	% correlation along y-axis
	Ry  = exp(-theta*abs(y));
	Sy  = real(ifft(Ry));
%	Ty  = sqrt(Sy);
	% why Sy and not sqrt Sy?
	Ty  = Sy;
	dx = x(2)-x(1);

	z = zeros(n,1);
	% columnwise
	for idx=1:n
	  for jdx=1:m
		% innovation
		e = zeros(n,1);
		if (idx==1 && jdx ==1)
			e(1) = 1;
		end
		%e = randn(m*n,1);
		% correlate along y
		e = real(ifft(Ty.*fft(e)));
		% advance the ornstein-uhlenbeck (AR1) process one column
		z = z + theta*(dx/m)*(mu_z-z) + sqrt(2*theta*dx/m)*sd_z*e;
		% state of the geometric process
	%	ai = exp(z);
		ai = z;
		% average m-rows
		ai = mean(reshape(ai,m,n),1)';
		% accumulate
		a(:,idx) = a(:,idx)+ai;
	end
	% normalize
	a = a./m^2;
end
end

