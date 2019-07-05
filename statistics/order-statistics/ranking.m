% So 2. Aug 11:56:41 CEST 2015
% Karl Kastner, Berlin
%
%% ranking for spearman statistics
%
function [ranking t] = ranking(x,w)
	if (isvector(x))
		x = cvec(x);
	end
	n = size(x,1);
	if (nargin()<2||isempty(w))
		w = (1:n)';
	else
		n
		w
		w = n*w/sum(w);
		%w = cumsum(w);
		w
	end
	m = size(x,2);
	[void ranking] = sort(x);
	w = w(ranking);
%	inverse = ranking(ranking)	
	lastid  = ones(1,m);
	t = zeros(1,m);
	r_ = ones(1,m);
	for jdx=1:m
	for idx=2:n
		if (void(idx,jdx) ~= void(lastid(jdx),jdx))
%			r_(lastid(jdx):idx-1,jdx) = mean(w(lastid(jdx):idx-1));
%			r_(lastid(jdx):idx-1,jdx) = sum(w(1:lastid(jdx)-1)) + mean(w(lastid(jdx):idx-1));
			r_(lastid(jdx):idx-1,jdx) = sum(w(1:lastid(jdx)-1)) + 0.5*sum(w(lastid(jdx):idx-1));
			if (idx-lastid(jdx) > 1)
				t(jdx) = t(jdx)+(idx-lastid(jdx));
			end % if
%			t = t+(idx-lastid)-1;		
			lastid(jdx) = idx;
		end % if
	end % idx
	% last row
%	r_(lastid(jdx):idx,jdx) = mean(w(lastid(jdx):idx));
	r_(lastid(jdx):idx,jdx) = sum(w(1:lastid(jdx)-1)) + 0.5*sum(w(lastid(jdx):idx));
	if (idx > lastid(jdx))
		t(jdx) = t(jdx)+(idx-lastid(jdx))+1;
	end % if
	%n = r_(end);
		ranking(ranking(:,jdx),jdx) = r_(:,jdx);
	end % jdx
	ranking
end

