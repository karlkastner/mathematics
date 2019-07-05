% Fri  2 Nov 13:30:02 CET 2018
% Karl Kastner, Berlin
%%
%% the product has two frequency components
%%
%%  product of the imaginary part of one and the real part of a second
%%  complex exponential
%%
%% input :
%% 	c : complex amplitudes
%%	o : frequencies
%% output :
%%	cp : amplitude of the product
%%	op : frequencies of the product
function [cp,op] = complex_exp_product_re_im(c,o)
	c(2)     = -1i*c(2);
	[cp,op]  = complex_exp_product_re_re(c,o);
end

