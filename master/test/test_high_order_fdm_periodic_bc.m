higher order bc with periodic bx: two wrapping grids : periodic (-1/+1); one test with extrapolated, three tests with extrapolated (+,-,simplectic); >> k=10; s=pi^2*(k:-1:1).^2'; n=100; A = (n+1)^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n); B=A; A(1,end) = -1; A(end,1) = -1; [Va Ea] = eigs(A,k,'SM'); [Vb Eb] = eigs(B,k,'SM'); subplot(3,1,1); plot(Va); subplot(3,1,2); plot(Vb); subplot(3,1,3); plot(Va.^2-Vb.^2) %)] + s*[1 1]

%
