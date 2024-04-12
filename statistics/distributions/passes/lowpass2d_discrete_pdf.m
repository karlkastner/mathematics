% Wed 27 Apr 11:02:07 CEST 2022
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
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% function [S,R,x,y] = lowpass2d_pdf_discrete_pdf(L,n,La,p);
function [S,R,x,y] = lowpass2d_pdf_discrete_pdf(L,n,La,p);
	% autocorrelation
	[R,x,y] = lowpass2d_discrete_acf(L,n,La);

	% spectal density
	S = fft2(R);

	% should be real up to rounding error, just to make sure
	S = real(S);

	if (nargin()>3 && ~isempty(p))
		S = S.^p;
	end
	% normalize to 1 at 0
	S = S/S(1,1);
end

