% Wed 27 Apr 11:05:35 CEST 2022
% Karl Kästner, Berlin
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%% function [S,R,r] = bandpass2d_discrete_pdf(L,n,Lf,p,q);
%
%% two dimensional spectral density of a discrete bandpass filter (dx finite)
function [S,R,r] = bandpass2d_discrete_pdf(L,n,Lf,p,q);
	if (nargin()<5)
		q = 1;
	end
	
	[S,R,r] = lowpass2d_discrete_pdf(L,n,Lf);
	sS = S.^(1/q);
	% bandpass density
	sS = sS.*(1-sS);
	S  = sS.^q;
	% higher order
	if (nargin()>3 && ~isempty(p))
		S = S.^p;
	end
	% normalize
	% note that the integral for p = 1 is 2 pi L^2
	df = 1./L;
	I = sum(S,'all')*df(1)*df(2);
	S = 2*S./I;
end

