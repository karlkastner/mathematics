>> f = expand(x*(x + y*(1 + 2*y)))
 
f =
 
x^2 + 2*x*y^2 + x*y
 
>> horner(f)
 
ans =
 
x*(2*y^2 + y + x)

