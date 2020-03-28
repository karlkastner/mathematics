n=100;
A=randn(n);
[Tv Tw D V W] = l2(A);
norm(A*V-V*Tv)
[T V W] = l3(A);
norm(V*T*W' - A)
