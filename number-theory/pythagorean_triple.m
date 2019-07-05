% 2015-08-25 12:18:03.498427818 +0200
% Karl Kastner, Berlin
%
%% pythagorean triple
%
function pt(n)
	x = (1:n)';
	XX = x*ones(1,n);
	
	XX2 = XX.^2 + XX'.^2;
	XX = sqrt(XX2);
	fdx = (0==mod(XX,1))
	[i j] = find(fdx);
	k = sqrt(i.*i+j.*j)
	[i j k k.*k] 
end
