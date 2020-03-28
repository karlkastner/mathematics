m=1000; N=10:10:100; for idx=1:length(N); n=N(idx); tic(); Test_jama.time_inverse(n, m); T(idx,1)=toc(), end; plot(N,T/m);

