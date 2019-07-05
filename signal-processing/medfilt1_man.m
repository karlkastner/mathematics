% Di 26. Mai 08:49:39 CEST 2015
% Karl Kastner, Berlin
%% moving median filter, supports columnwise operation
function [Y S L U] = medfilt1_man(X,n)
	Y = zeros(size(X));
	S = zeros(size(X));
	L = zeros(size(X));
	U = zeros(size(X));
	l = floor(n/2);
	r = round(n/2)-1;
	for idx=1:size(X,1)
		[Y(idx,:) S(idx,:) L(idx,:) U(idx,:)] = median_man( ...
			X(max(1,idx-l):min(end,idx+r),:)');
	end
end % medfilt1_man

