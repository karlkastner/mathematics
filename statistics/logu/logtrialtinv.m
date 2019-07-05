% Fri 23 Mar 13:35:25 CET 2018
%
%% inverse of the logarithmic triangular distribution
%
function x = logtriinv(a,b,c,F)
	e = exp(1);
	d    =     (b*(log(c) - log(b) + 1))/((log(b) - log(c))) ...
                 - (a - b - b*log(a) + b*log(b))/((log(a) - log(b))) ...
                 - c/((log(b) - log(c)));
%	F   = F*d;
%	whos
%	x   = (-a-F)./lambertw(0,(-a-F)/(a*e));
	x = (d*F*log(b/a)-a)./lambertw_numeric((-a - d*F*log(a/b))/(a*e),0);

	Fb  = -(a - b - b*log(a) + b*log(b))/(d*(log(a) - log(b)));
	fdx = F > Fb;
	% TEST: logtricdf(a,b,c,b)-Fb

%	sum(fdx)/length(fdx)
%	F_ = F(fdx);
%	x(~fdx) = NaN;
%if (1)
%	x(fdx) = (F_*log(a)*log(b) + a*log(b) - b*log(a) - F_*log(a)*log(c) - a*log(c) + F_*log(b)*log(c) + b*log(c) - F_*log(b)^2) ...
%		./ ((log(a) - log(b)) ...
%		    * lambertw(0,(a^(-1/(log(a) - log(b)))*(b^(-log(c)/log(a) - 1/log(a))*c)^(-log(a)/(log(a) - log(b))) ...
 %                                               *(-F_*log(b)^2 + a*log(b) + F_*log(a)*log(b) + F_*log(c)*log(b) - b*log(a) - a*log(c) + b*log(c) - F_*log(a)*log(c))) ...
  %                     / (log(a) - log(b))));
%	F(fdx) = Fc

	%F = F/d;
	F_ = F(fdx) - (b*(log(c) - log(b) + 1))/(d*(log(b) - log(c))) ...
		    + (a - b - b*log(a) + b*log(b))/(d*(log(a) - log(b)));
	
	x(fdx) =  (d*F_*log(b/c))./lambertw_numeric((d*F_*(log(b/c)))/(c*e),-1);

%if (0)
%	x(fdx) = (d*F(fdx)*log(a)*log(b) + a*log(b) - b*log(a) - d*F(fdx)*log(a)*log(c) - a*log(c) + d*F(fdx)*log(b)*log(c) + b*log(c) - d*F(fdx)*log(b)^2) ... 
%			./ ((log(a) - log(b))*lambertw(0,(a^(-1/(log(a) - log(b)))*(b^(-log(c)/log(a) - 1/log(a))*c)^(-log(a)/(log(a) - log(b))) ...
%                                                         .*(-d*F(fdx)*log(b)^2 + a*log(b) + d*F(fdx)*log(a)*log(b) + d*F(fdx)*log(c)*log(b) - b*log(a) - a*log(c) + b*log(c) - d*F(fdx)*log(a)*log(c))) ...
%                                                         /(log(a) - log(b))));
%end
%	x
%% = (d F log(a) log(b) + a log(b) - b log(a) - d F log(a) log(c) - a log(c) + d F log(b) log(c) + b log(c) - d F log^2(b))/((log(a) - log(b)) W((a^(-1/(log(a) - log(b))) (b^(-log(c)/log(a) - 1/log(a)) c)^(-log(a)/(log(a) - log(b))) (-d F log^2(b) + a log(b) + d F log(a) log(b) + d F log(c) log(b) - b log(a) - a log(c) + b log(c) - d F log(a) log(c)))/(log(a) - log(b))))
%%x = (d F log(a) log(b) + a log(b) - b log(a) - d F log(a) log(c) - a log(c) + d F log(b) log(c) + b log(c) - d F log^2(b))/((log(a) - log(b)) W((a^(-1/(log(a) - log(b))) (b^(-log(c)/log(a) - 1/log(a)) c)^(-log(a)/(log(a) - log(b))) (-d F log^2(b) + a log(b) + d F log(a) log(b) + d F log(c) log(b) - b log(a) - a log(c) + b log(c) - d F log(a) log(c)))/(log(a) - log(b))))
%
%
%
%else
%	x(fdx) = NaN;
%end
	% todo check true f
	x(F<0) = NaN;
	x(F>1) = NaN;
end

