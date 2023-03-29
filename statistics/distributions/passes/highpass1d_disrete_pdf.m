% Sat 26 Jun 21:04:27 CEST 2021
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
%
%% spectral density of the pth-order high-pass
%%
%% Note that there are two alternative definitions
%% S_hp = S_h1^p = (1 - S_l1)^p (recursive highpass-filtering)
%% or S_hp~ = (1 - S_lp^p) (1 - recursive lowpass-filtering)
%% here, recursive highpass filtering is represented
%%
%% Sh = |F*Ah|^2p
%% Sh^(1/2p) = F Ah = F(I-Al)
%%           = F(I - Al)
%%           = F(I - F^-1 Sl^1/2 F)
%%           = (F F^-1 - Sl^1/2) F
%%           = (I - Sl^(1/2)) F 
%%
%% function [S_h,S_h1] = spectral_density_hp(f,r0,fmax,order,varargin)
function [S_hp, S_h1] = highpass1d_discrete_pdf(fx,arg2,order,dx,varargin)
	if (nargin()<3||isempty(order))
		order = 1;
	end
	if (length(varargin)>1)
		normalize = varargin{2};
		varargin  = varargin(1);
	else
		normalize = true;
	end
	if (length(varargin)<1)
		form = 'rho';
	else
		form = varargin{1};
	end

	%r      = highpass_arg(arg2,dx,order,form);
	%r      = highpass_arg(arg2,dx,1,form);
	r      = bandpass_arg(arg2,dx,varargin{:});
	[S_lp,S_lp1]   = lowpass1d_discrete_pdf(fx,r,1,dx,'r',false);
	S_hp     = (1-S_lp1.^0.5).^(2*order);
	% by default, normalize at least to 1 for fmax
	if (normalize > 0)
		[S_lp_max,S_lp1_max] = lowpass1d_discrete_pdf(0.5/dx,r,1,dx,'r',false);
		S_hp_max = (1-S_lp1_max.^0.5).^(2*order);
		S_hp = S_hp/S_hp_max; 
	end
	if (~issym(fx) && (normalize>0))
		% normalize
		%df = 1/L;
		%S_bp = 2*S_bp/sum(S_bp*df);
		I = spectral_density_area(fx,S_hp);
		S_hp = S_hp./I;
	end
end
 
