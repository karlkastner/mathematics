% Mon  9 Jan 11:12:08 CET 2017
% Karl Kastner, Berlin
%
%% weighted mean of angles
%
% input
% alpha : x*m, [rad] angle
%
% ouput
% ma    : 1*m, [rad] mean angle
% sa    : 1*m, [rad] standard error of mean angle for uncorrelated error
function [ma sa] = meanangle(alpha,dim,w)
	if (nargin()<2)
		dim = 1;
		if (isvector(alpha))
			alpha = cvec(alpha);
		end
	end
	n = size(alpha,dim);

	s = sin(alpha);
	c = cos(alpha);

	if (nargin()<3)
		ms = nanmean(s,dim);
		mc = nanmean(c,dim);
	else
		ms = wmean(w,s,dim);
		mc = wmean(w,c,dim);
	end

	% mean angle
	ma = atan2(ms,mc);
	
	% residual
	alpha = wrapToPi(alpha);
	res = bsxfun(@minus,alpha,ma);
	res = wrapToPi(res);
	
%	squeeze(nanrms(bsxfun(@minus,s,ms),dim))
%	squeeze(nanrms(bsxfun(@minus,c,mc),dim))
%	squeeze(nanrms(res,dim))

	% standard error
	% sa = serr(res);
	% TODO weighing
	% TODO n needs to be adapted to number of valid samples per column
	sa = sqrt(nanmean(res.^2,dim)/((n-1)));
end % mean_angle

