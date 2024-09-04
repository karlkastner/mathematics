% 2024-06-30 11:54:55.114276491 +0200
function [f0,s] = cauchymirroredpdf_mode2par(fc,Sc,par0)
	if (fc == 0)
		f0 = NaN;
		s  = NaN;
	else
	f0 = (fc*((16*Sc^2*fc^2*pi^2 + 1)*(16*Sc^2*fc^2*pi^2 - 3))^(1/2))/(16*Sc^2*fc^2*pi^2 - 1);
	s = (8*Sc*fc^2*pi)/(16*Sc^2*fc^2*pi^2 - 1);	 % 1 3
%	s = -((-(2*pi*Sc*fc - 1)*(2*pi*Sc*fc + 1))^(1/2) - 1)/(2*Sc*pi);
%	s = ((-(2*pi*Sc*fc - 1)*(2*pi*Sc*fc + 1))^(1/2) + 1)/(2*Sc*pi);
%	s = (8*Sc*fc^2*pi)/(16*Sc^2*fc^2*pi^2 - 1);	 % 1 4
%	s = -((-(2*pi*Sc*fc - 1)*(2*pi*Sc*fc + 1))^(1/2) - 1)/(2*Sc*pi); % 3 4
%	s = ((-(2*pi*Sc*fc - 1)*(2*pi*Sc*fc + 1))^(1/2) + 1)/(2*Sc*pi); % 3 4


%	   r = 2*fc*sqrt(fc^2 + s^2) - fc^2 - s^2;
%	   if (r >= 0)
%	   	f0 = sqrt(r);
%	   else
%	   	f0 = sqrt(- 2*fc*sqrt(fc^2 + s^2) - fc^2 - s^2);
%	   end
%	 f0 = (fc^2/3 - (2*((fc^2 + fc*s + s^2)*(fc^2 - fc*s + s^2))^(1/2))/3 - s^2/3)^(1/2)
% (fc^2/3 - (2*((fc^2 + fc*s + s^2)*(fc^2 - fc*s + s^2))^(1/2))/3 - s^2/3)^(1/2)
%	f0 = sqrt((fc^2 - 2*sqrt((fc^2 + fc*s + s^2)*(fc^2 - fc*s + s^2)) - s^2)/3);
%	f0 = ((s + s*((- 16*s*Sc^2*fc^2*pi^2 + 8*Sc*fc^2*pi + s)/s)^(1/2) + 2*pi*Sc*fc^2 - 2*pi*Sc*s^2)/(2*Sc*pi))^(1/2);

%f0= ((fc*1i + (fc*(- Sc^2*fc^3*pi^2*16i + 8*Sc*fc^2*pi + fc*1i)*1i)^(1/2) + 4*pi*Sc*fc^2)/(2*Sc*pi))^(1/2);
%                                                                                           (fc^2)^(1/2)
%f0=-((fc*1i + (fc*(- Sc^2*fc^3*pi^2*16i + 8*Sc*fc^2*pi + fc*1i)*1i)^(1/2) + 4*pi*Sc*fc^2)/(2*Sc*pi))^(1/2)
%                                                                                         -(fc^2)^(1/2)
%f0= (-(fc*1i + (-fc*(Sc^2*fc^3*pi^2*16i + 8*Sc*fc^2*pi - fc*1i)*1i)^(1/2) - 4*pi*Sc*fc^2)/(2*Sc*pi))^(1/2)
%f0=-(-(fc*1i + (-fc*(Sc^2*fc^3*pi^2*16i + 8*Sc*fc^2*pi - fc*1i)*1i)^(1/2) - 4*pi*Sc*fc^2)/(2*Sc*pi))^(1/2)


%	 ((s - s*((- 16*s*Sc^2*fc^2*pi^2 + 8*Sc*fc^2*pi + s)/s)^(1/2) + 2*pi*Sc*fc^2 - 2*pi*Sc*s^2)/(2*Sc*pi))^(1/2)

%	Sc = cauchymirroredpdf(fc,f0,
%	sqrt(f0^2 - fc^2);

% ((2*((f0^2 + f0*s + s^2)*(f0^2 - f0*s + s^2))^(1/2))/3 + f0^2/3 - s^2/3)^(1/2)
%-(f0^2/3 - (2*((f0^2 + f0*s + s^2)*(f0^2 - f0*s + s^2))^(1/2))/3 - s^2/3)^(1/2)
%-((2*((f0^2 + f0*s + s^2)*(f0^2 - f0*s + s^2))^(1/2))/3 + f0^2/3 - s^2/3)^(1/2)
	end

if (0)
	if (nargin()<3)
		[par(1),par(2)] = cauchypdf_mode2par(fc,Sc);
	end
	par = lsqnonlin(@(par) [fc,Sc] - fun(par),par);
	f0 = par(1);
	s  = par(2);
end
	function out = fun(in)
		[fc,Sc] = cauchymirroredpdf_mode(in(1),in(2));
		out = [fc,Sc];
	end
end
