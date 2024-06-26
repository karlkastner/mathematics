% 2018-11-02 13:42:52.809075592 +0100
% Karl Kastner, Berlin

o = rand(2,1);
o = [1,2];
c  = rand(2,1) + 1i*rand(2,1);
t = linspace(0,2*pi)';

id = 2;
switch (id)
case {0}
[cp,op] = complex_exp_product_re_re(c,o);

y = [real( c(1)*exp(1i*o(1)*t) ).*real(c(2)*exp(1i*o(2)*t)) ...
	(real(cp(1)*exp(1i*op(1)*t) ) + real(cp(2)*exp(1i*op(2)*t)) ) ];
case {1}
[cp,op] = complex_exp_product_re_im(c,o);
y = [real( c(1)*exp(1i*o(1)*t) ).*imag(c(2)*exp(1i*o(2)*t)) ...
	(real(cp(1)*exp(1i*op(1)*t) ) + real(cp(2)*exp(1i*op(2)*t)) ) ];
case {2}
[cp,op] = complex_exp_product_im_im(c,o);
y = [imag( c(1)*exp(1i*o(1)*t) ).*imag(c(2)*exp(1i*o(2)*t)) ...
	(real(cp(1)*exp(1i*op(1)*t) ) + real(cp(2)*exp(1i*op(2)*t)) ) ];
end

figure(1)
clf
plot(t,y)

norm(y(:,1) - y(:,2))

syms c_1 c_2 o_1 o_2
[cp,op] = complex_exp_product_re_re([c_1,c_2],[o_1,o_2])
