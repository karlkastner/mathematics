syms a b c d e f g h i l
	A = reshape([a,b,c,d,e,f,g,h,i],3,3)
ev = eig(A)
A = reshape([a-l,b,c,d,e-l,f,g,h,i-l],3,3)
 collect(det(A),'l')



