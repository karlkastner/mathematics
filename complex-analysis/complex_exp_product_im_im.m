% Fri  2 Nov 13:30:02 CET 2018
% Karl Kastner, Berlin
%
%%  product of the imaginary part of two complex exponentials
%%
%% the product has two frequency components
%%
%% input :
%% 	c : complex amplitudes
%%	o : frequencies
%% output :
%%	cp : amplitude of the product
%%	op : frequencies of the product
function [cp,op] = complex_exp_product_im_im(c,o)
	c        = -1i*c;
	[cp,op]  = complex_exp_product_re_re(c,o);
end

