function derive_sum

% derive sum formula

% sum i
% s = a2 n^2 + a1 n + a0
A = [ 1 1 1
      4 2 1
      9 3 1];
% cumsum
F = [ 1 3 6 ]';
s=A \ F
m = 1/s(1)
s = m*s 

% sum sum sum i
% s = a4 n^4 + a3 n^3 + a2 n^2 + a1 n + a0
A = [   1   1  1  1  1
       16   8  4  2  1
       81  27  9  3  1
      256  64 16  4  1
      625 125 25  5  1 ]
F = [ 1     5    15    35    70 ]';
s = A \ F
m = 1/s(1)
s = m*s 
%s = round(s*1e12)*1e-12;
%[n d] = numden(sym(s))

n=10;
L=1;
X = sum1(1:n/2+1) - 0.5
X = L*[-fliplr(X) X]/X(end)
diff(X)
plot(X,ones(length(X)),'.')
%length(find(abs(X)<1/n))
%length(find(abs(X)>1/n))
%harmmean(sum1(1:n)/sum1(n))
%harmmean(sum3(1:n)/sum3(n))
1/sum3(n)
end % derive_sum

