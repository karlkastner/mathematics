% Thu  4 May 17:29:38 CEST 2017
% Karl Kastner, Berlin
%
%% nth-roots of a complex number
%%
%% input:
%% c : complex number
%% n : order of root
%%     n must be rational, to obtain n solutions
%%     otherwise no finite set of solutions exists
%%
%% r : roots of the complex number
function r = croots(c,n)
	if (~issym(c))
		r = roots([1, zeros(1,n-1), -c]);
	else
		% transform roots to power
		r = (abs(c))^(1/n)*exp(1i*angle(c)/n + 1/n*2i*pi*(0:n-1)');
	end
end

