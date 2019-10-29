% Fri Nov  4 18:46:23 MSK 2011
% Karl Kästner, Berlin

function D1 = gradient(X)
	n = length(X);
	xl = X(1:end-1);
	xr = X(2:end);
	D1 = spdiags(0.5./(xr - xl)*[-1 0 1],-1:1,n-2,n-2);
end % function gradient

% \begin{align}
% u(y)'' + (V(y) - l) u(y) = 0
% y = f(x)
% u(y)'' + (V(y) - l) u(y) = 0
% ( y'u'(y) )' + (V(y) - l) u(y) = 0
%  y''u'(y) + (y')^2 u''(y)  + (V(y) - l) u(y) = 0
% (diag(D_x^2 F x)D_u + diag(D_x F x)^2 D_u^2) u(y) = l u(y)
% D_x : different bc : [x1-h | x1 .. xn | xn+h ]
% D_u : original constant step width or new variable step width? assume the letter
% \end{align}
% test l and v for V=0 and f(x) = x^2
% => 2 I D_u + 4 diag(x^2) D_u^2 ~ D_u^2 in original grid

