% Fri  2 Nov 13:30:02 CET 2018
% Karl Kastner, Berlin
%
%%  product of the imaginary part of one and the real part of a second
%%  complex exponential
%%
%% the product has two frequency components
%%
%% input :
%% 	c : complex amplitudes
%%	o : frequencies
%% output :
%%	cp : amplitude of the product
%%	op : frequencies of the product
function [cp,op] = complex_exp_product_im_re(c,o)
	c(1)     = -1i*c(1);
	[cp,op]  = complex_exp_product_re_re(c,o);
end

