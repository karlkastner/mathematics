% Di 2. Feb 10:17:38 CET 2016
% Karl Kastner, Berlin
%% triangular filter window
function win = triwin(nf)
	nfh = ceil(nf/2);
	win = [1:nfh,nf-nfh:-1:1]';
	win = win/sum(win);
end

