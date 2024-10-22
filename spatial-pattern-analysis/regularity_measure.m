% 2024-07-01 13:09:10.093941004 +0200
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
%% compute various regularity measures from the spectral density,
%% all retransformed to the scale-invariant measure Sc/lambda_c
%
function [reg, leg_C, out] = regularity_measure(fx,S,distribution)
	fdx = (fx>=0);
	% spectral resolution
	df  = fx(2)-fx(1);
	% spatial extent
	L = 1./df;
	% number of grid points
	n = length(fx);
	% spatial axis
	x = (0:n-1)'*L/n;

	% normalize
	S = S./sum(S.*df);
	% density along positive half-plane
	Sp = 2*S.*(fx>=0);

	% autocorrelation function (periodically extended)
	R  = real(fft(S));
	R  = R/R(1);

	% mode of S
	[Spc,mdx] = max(Sp);
	fc = fx(mdx);
	if (1)
		[Spc, fc] = extreme3(fx,Sp,mdx);
	end

	% cumulative density
	C    = cvec(cumsum(Sp(fdx)))*df;
	fx_  = fx(fdx);
	fdx_ = [true; C(2:end) ~= C(1:end-1)];
	% median frequency
	fme = interp1(C(fdx_),fx_(fdx_)+0.5*df,0.5,'linear');
	% quartiles 
	f25 = interp1(C(fdx_),fx_(fdx_)+0.5*df,0.25,'linear');
	f75 = interp1(C(fdx_),fx_(fdx_)+0.5*df,0.75,'linear');

	%fc_ = fc;
	fc_ = fme;

	% moments of the density
	fmu = sum(Sp.*fx)*df;
	fs2 = sum(Sp.*(fx-fmu).^2)*df;
	fsd = sqrt(fs2);

	% right point where the density drops to half its max value
	idr = find( (Sp <= 0.5*Spc) & (fx>fc) & (fdx), 1, 'first');

	% interpolate linearly
	% S = S(idr-1) + (f-fx(idx-1))*(S(idr)-S(idr-1))/df;	
	% fr  = fx(idr);
	if (isempty(idr))
		fr = max(fx);
	else
		fr =  (0.5*Spc - Sp(idr-1))*df/(Sp(idr)-Sp(idr-1)) + fx(idr-1);
	end

	% left point where the density drops below half its maximum value
	idl = find( (Sp <= 0.5*Spc) & (fx<fc) & (fx>=0), 1, 'last');

	if (isempty(idl))
		fl = 0;
	else
		% S = S(idl) + (f-fx(idl))*(S(idl+1)-S(idl))/df;
		fl = (0.5*Spc - Sp(idl))*df/(Sp(idl+1)-Sp(idl)) + fx(idl);
		%fl  = fx(idl);
	end
	% width of the distribution
	w = fr - fl;

	% entropy
	lnS = log(Sp);
	lnS(Sp==0) = 0;
	h = -sum(lnS(fdx).*Sp(fdx))*df;

	% mode of acf
	xdx = (x>0.5*fc & x < 1.5*fc);
	[Rc, rcdx]   = max(R.*xdx);
	if (1)
		[Rc, xc] = extreme3(x,R,rcdx);
	end

	%fc = fme;

	% mode Sc, reg = Sc/lc
	reg(1) = Spc.*fc;

switch (distribution)
case {'gauss','normal'}
% standard deviation s
%	Sc = 1/(sqrt(2 pi) sd)
%	reg = 1/(sqrt(2*pi)*sd*lc)
	reg(2) = fc./(sqrt(2*pi)*fsd);
% width
	reg(3) = fc./(1./sqrt(log(256))*sqrt(2*pi)*w);
	% interquartile range
	reg(4) = (erfinv(0.5)-erfinv(-0.5)).*fc_./(sqrt(pi)*(f75-f25));
% entropy h   = ln(sqrt(2 pi)s) + 1/2
%	  s   = exp(h - 1/2)/sqrt(2 pi)
%	  Sc  = 1/exp(h - 1/2)
%	  reg = 1/exp(h - 1/2)*lc)
	  reg(5) = fc.*exp(0.5-h);
	 % reg(idx,3) = fc./exp(-(h-0.5));
% maximum acf
%	  Rc  ~ exp(-2*pi^2*s^2*lc^2)
%	  s*lc = sqrt(-log(Rc)/2)/pi
%	  reg  = sqrt(-pi/(log(Rc)))
	  reg(6) = sqrt(pi/log(1/Rc));
case {'laplace'}
	% mode Sc = 1/(2s)
	% standard deviation s
	% sd = sqrt(2)s
	% s = sd/sqrt(2)
	% Sc/lc = 1/(sqrt(2)sd*lc)
	reg(2) = fc./(sqrt(2)*fsd);
% width % TODO
	% w =  -2 ln(1/2) s
	% Sc = 1/(-1/ln(1/2) w)
	reg(3) = fc*log(2)/w;
	% interquartile range
	reg(4) = log(2).*fc_/((f75 - f25));
% entropy
	% h = ln (2s) + 1
	% exp(h - 1)/2
	% reg = 1./(exp(h-1)*lc)
	reg(5) = fc.*exp(1-h);
% mode of acf
	% s lc = sqrt((1/Rc - 1)/(4 pi^2))
	reg(6) = pi.*sqrt(Rc./(1 - Rc));
case {'cauchy','lorentz'}
	% std
	reg(2) = NaN;
	% width
	% w = 2*s
	% s = w/2
	reg(3) = 2/(pi*w)*fc;
	% interquartile range
	reg(4) = (tan(pi/4)-tan(-pi/4))*fc_./(pi*(f75 - f25));
	% entropy
	% h = ln(4 pi s)
	% s = exp(h)/(4 pi)
	reg(5) = 4*exp(-h)*fc;
% mode of acf
	%Rc = exp(-2 pi s lc)
	%s*lc = -log(Rc)/(2*pi)
	% Sc  = 1/(pi*-log(Rc)/(2*pi))
	reg(6) = 2/log(1/Rc);
end
	out.Spc      = Spc;
	out.fc      = fc;
	out.Rc      = Rc;
	out.entropy = h;
	out.mu      = fmu;
	out.fsd      = fsd;
	out.width   = w;
	out.f25     = f25;
	out.f75     = f75;
	out.fme     = fme;

	leg_C = {'Mode','Standard dev.','Width','IQR','Entropy','Autocorrelation'};
end % regularity measure

