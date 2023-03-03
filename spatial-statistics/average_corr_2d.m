% Tue 28 Feb 14:45:20 CET 2023
% requires equispaced grid
function acor = average_corr_2d(fun,x,y,order)
	dx = x(2)-x(1);
	dy = y(2)-y(1);
	xx = flat(repmat(cvec(x),1,length(y)));
	yy = flat(repmat(rvec(y),length(x),1));
	[w,p] = int_1d_gauss(order);
	% integrate along x for x1
	x1 = 0   + p*(0.5*dx*[-1;1])
	% integrate along x for x2;
	x2 = xx + (p*(0.5*dx*[-1;1]))';
	% integrate along x for y1
	y1 = 0   + p*(0.5*dy*[-1;1])
	% integrate along x for y2
	y2 = yy + (p*(0.5*dy*[-1;1]))';
	acor = zeros(numel(xx),1);
	for x1id=1:length(w)
	 for x2id=1:length(w)
          for y1id=1:length(w)
           for y2id=1:length(w)
		% note we do not need to multiply and divide by dx and dy,
		% as they multiplication for integration cancels with division for averaging
		acor = acor + w(x1id)*w(x2id)*w(y1id)*w(y2id)*fun(x1(x1id),x2(:,x2id),y1(y1id),y2(:,y2id));
	   end
	  end
	 end
	end
	acor = reshape(acor,length(x),length(y)); 
end % average_corr_2d

