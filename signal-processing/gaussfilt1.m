% Mi 8. Apr 18:19:56 CEST 2015
% Karl Kastner, Berlin
%% filter data series with a gaussian window
function [X f] = gaussfilt1(X,n)
%	f = ones(n,1)/n;
	m = 7;
	f = normpdf((-m*n:m*n)/(2*n))';
	f = f/norm(f);
	
	X = wmeanfilt(f,X);
	%X = winfilt1(f,X);
%	if (isvector(X))
%		%X = conv(X,f,'same');
%		X = conv1_man(cvec(f),cvec(X));
%	else
%		for idx=1:size(X,2);
%			X(:,idx) = conv(X(:,idx),f,'same');
%		end
%	end
end

