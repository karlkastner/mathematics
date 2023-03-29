% Fri  3 Mar 13:08:40 CET 2023
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
function [fc,Sc] = bandpass2d_pdf_mode(f0,order,L,n)
	if (1)
		m = 100;
		L = m./(f0/8);
		n = m^2;
		%S = bandpass2d_pdf_hankel(fx,par(1),par(2));
		df = 1./L;
		fx = fourier_axis(L,n);
		%S  = bandpass2d_pdf_exact(abs(fx),f0,order);
		S = bandpass2d_pdf_hankel(L,n,f0,order);
		% normalize
		%S = 2*S./sum(S*df);
		%S_ = 2*S_./sum(S_*df);
		[Sc,mdx] = max(S);
		if (0)
			fc = fx(mdx);
		%fc = abs(fx(mdx));
		else
			[Sc,fc] = extreme3(fx,S,mdx);
		end
%		par = par
	end

	if (0)
		%nargout()>1)
	%	fc0 = 1.0;
	%	fc = lsqnonlin(@(f) (bandpass2d_pdf(f,f0,order)-1.0),fc0);
	
		n  = 20;
		l0 = 6/f0;
		L  = n*l0;
		df = 1/L;
		fx = fourier_axis(L,n*n)';
		S = bandpass2d_pdf(abs(fx),f0,order);
		S = 2*S./sum(S*df);
		% TODO quadratic max
		[Sc,mdx] = max(S);
		fc = abs(fx(mdx));
		%mode = [fc,Sc]
	end
end

