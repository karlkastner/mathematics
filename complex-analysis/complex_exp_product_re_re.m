% Tue  1 May 10:13:21 CEST 2018
% Karl Kastner, Berlin
%%
%%  product of the real part of two complex exponentials
%%
%% re(c1 exp(io1x))*re(c2 exp(io2x)) = 
%%	1/2*(     real(c1*c2*exp(i*(n1+n2)*o*x)) ...                            
%%               + real(conj(c1)*c2*exp(i*(n2-n1)*o*x)) )
%%
%% the product has two frequency components
%%
%% input :
%% 	c : complex amplitudes
%%	o : frequencies
%% output :
%%	cp : amplitude of the product
%%	op : frequencies of the product
function [cp, op] = complex_exp_product_re_re(c,o)
	if (isnumeric(o) && o(1) > o(2))
		c = [c(2),c(1)];
		o = [o(2),o(1)];
	end
	cp = 1/2*[c(1)*c(2), conj(c(1))*c(2)];
	op = [o(1)+o(2), o(2)-o(1)];
end

