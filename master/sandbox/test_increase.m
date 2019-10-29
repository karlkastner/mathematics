	javaaddpath('.')
	A = reshape(1:6,2,3)
%	A = javaArray('double', 4, 5)
%	java.lang.double(10)
	Test.display(A);
	A = Test.increase(A);
	A

