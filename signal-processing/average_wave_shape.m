% Tue  5 Oct 21:56:26 CEST 2021
%
%% extract waves with varying length from a wave train and and average their shape
%
function [y_mu,y_std,yy,ps,pi_] = average_wave_shape(y,li)
	%method = 'linear';
	method = 'spline';
	single = false;
	y = y-min(y);
	y = y/max(y);
	if (single)
		thresh = graythresh(y)
	else
		mthresh = multithresh(y,2);
	end
	n = length(y);
	ps = [];
	mindx = 2;
	state = 1;
	for idx=2:n
if (single)
		if (y(idx)<y(mindx))
			mindx =idx;
		end
		if (y(idx-1)<=thresh && y(idx)> thresh)
			%ps(end+1) = mindx;
			ps(end+1) = idx;
			mindx = idx;
		end
else
		switch (state)
		case {1} % search for end
			if (y(idx)<mthresh(1))
				state = 2;
				mindx = idx;
			end
		case {2} % search for start
			if (y(idx)<y(mindx))
				mindx = idx;
			end
			if (y(idx)>mthresh(2))
				ps(end+1,1) = idx; 
				%ps(end+1,1) = mindx; %idx;
				state = 1;
			end
		otherwise
		end
end
	end
	% get average length
	d = 1;
	ll = ps(d+1:end)-ps(1:end-d);
	if (nargin()<2)
	li = round(2*mean(ll));
	end
	yy = zeros(li+0,length(ps)-d);
	pi_ = ((0:li-1)+0.5)/li;
	%w = [];
	for idx=1:length(ps)-d
		id = ((ps(idx)-1:ps(idx+d)+1));
		yy(:,idx) = interp1((id-ps(idx))/(ps(idx+1)-ps(idx)),y(id),pi_,method,'extrap');
		%yy(:,idx) = interp1((id-id(1))/ll(idx),y(id),pi_,method,'extrap');
		%w(idx,1) = sum(y(id(2:end-1)).^2);
	end
	% weights
	w = ll;
	%w = ones(length(ps)-1,);
	w = w./sum(w);
	%w(end+1) = 0;
	%w = w/sum(w);
%	tail      = y(mindx:end);
%	yy(1:length(tail),idx) = tail;
%	cc(1:length(tail),idx) = 1;
	y_mu  = yy*cvec(w);
	y_std = zeros(size(y_mu));
	for idx=1:size(yy,1)
		y_std(idx,1) = wstd(w,yy(idx,:));
	end
%	y_std = zeros(size(y_mu));
%	y_me = zeros(size(y_mu));
%	for idx=1:size(yy,1)
%		y_std(idx,1) = std(yy(idx,cc(idx,:)==1),[],2);
%		y_me(idx,1)= median(yy(idx,cc(idx,:)==1),2);
%	end
%	y_se = y_std./sqrt(sum(cc,2));
end
