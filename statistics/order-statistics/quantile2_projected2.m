% 2016-03-01 11:45:59.744766302 +0100
%% spatial qunatile for chosen direction
function Q = quantile2_directional(X,p,dir)

	
	me = median(X,2);
	dX = bsxfun(@minus,X,me);
	% scale
	scale = 0.5*diff(quantile(X',normcdf([-1 1]))',[],2);
	dX = bsxfun(@times,dX,1./scale);

	% if no direction specified, take X as diretion
	if (nargin() < 3)
		dir = dX;
	end
	% TODO avoid div by zero
	dir = bsxfun(@times,dir,1./sqrt(sum(dir.^2)));

	% for each direction
	for idx=1:size(dir,2)
		% project to angle i (first element of rotation matrix)
		if (0)
			R = [dir(1,idx),-dir(2,idx),
                	     dir(2,idx),dir(1,idx)];
			Xp = R'*dX;
			Xp=Xp(1,:);
		else
			Xp = dir(:,idx)'*dX;
			% projected quantile
			q = quantile(Xp,p);
		end
		% nearest point of projected quantile
		[void pdx] = min( (Xp-q).^2 );
		q_ = Xp(pdx);
		Q(:,idx) = q_*dir(:,idx);
	end	
	% scale back
	Q = bsxfun(@times,Q,scale);
	% translate back
	Q = bsxfun(@plus,Q,me);
end

