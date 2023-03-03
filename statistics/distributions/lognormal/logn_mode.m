% Thu  9 Sep 12:57:04 CEST 2021
% Karl Kästner, Berlin
%
%% mode (maximum) of the log-normal density
%
% function [x_mode,y_mode] = logn_mode(lmu,lsd)
function [x_mode,y_mode] = logn_mode(lmu,lsd)
	x_mode = exp(lmu - lsd.^2);
	[mu,sd] = logn_param2moment(lmu,lsd);
	y_mode = (1 + (sd./mu).^2)./(mu.*sqrt(2*pi*log(1 + (sd./mu).^2)));
end

