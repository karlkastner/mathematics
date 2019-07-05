% 2016-03-03 10:22:54.035514656 +0100
% Karl Kastner, Berlin
%% hanning filter window
function f = hanwin_(x,mu,L)
	dx = x - mu;
	f = 2./L.*cos(pi*dx./L).^2.*(dx > -L/2 & dx < L/2);
end

