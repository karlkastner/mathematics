% 2021-06-21 15:53:50.278752860 +0200
%
% radial "periodogram", this is actually a density estimate, for fr that are not too small
%
%function [S,fri,se,count] = periodogram_radial(Shat,L)
%
% Shat : 2-dimensional periodogram
% L    : domain length
% 
% S.mu : Radially averaged periodogram
%       S_mu(k)  = 1/k int_0^{2 pi} Shat(k,theta) d theta
% S.normalized :  normalized radially averaged periodogram,
%       i.e. the radial density estimate
%       Sn  = S/int_0^inf S df
% rS  : k*S = int Shat dtheta = k S
%
% when S is flattened into a vector, the isotropic part of the 2D density can be recovered with:
% S_iso     = (A*S_radial)
% S_radial  = A^-1 S_hat
%
function [S, fri, count, A] = periodogram_radial(Shat,L)
	n = size(Shat);
	if (nargin()<2)
		L = n;
	end
	fx = fourier_axis(L(1),n(1));
	fy = fourier_axis(L(2),n(2));
	fmax = max(max(fx),max(fy)); 
	df   = 1./L;
	dfi  = sqrt(df(1).*df(2));
	fr   = hypot(fx,fy');
	frmax = hypot(max(fx),max(fy));
	fri = (0:dfi:ceil(frmax+dfi))';
%	fri  = (0: 

	% this is a linear interpolation
	% each value is proportionally split between the next lower and next larger bin
	% integer part
	rat = fr/dfi;
	% index into 1D radial density
	id_Sr = floor(rat);
	p   = 1-(rat-id_Sr);
	% 0 is id_Sr 1
	id_Sr = id_Sr+1;

	s = [length(fri)+1,1];
	S  = struct();

	% half-sum
	sumS = 0.5*( accumarray(id_Sr(:),p(:).*Shat(:),s,@sum) ...
	         + accumarray(id_Sr(:)+1,(1-p(:)).*Shat(:), s, @sum) );
	% number of id_Srs
	count = ( accumarray(id_Sr(:),p(:),s,@sum) ...
	        + accumarray(id_Sr(:)+1,(1-p(:)), s, @sum) );
	sumS  = sumS(1:end-1);
	count = count(1:end-1);

	%n  = accumarray(flat(r),ones(numel(f),1),[],@sum);
	S.mu = sumS./count;
	S.mu(count==0) = 0;
	% normalized
	S.normalized = S.mu./(sum(mid(S.mu)).*dfi);

	% matrix
	if (nargout()>3)
		nr = length(fri);
		nn = numel(Shat);
		% index into 2D periodogram
		id_S     = (1:nn)';
%		id_S = reshape(id_S,size(Shat));
%		id_S = id_S';
%		id_S = id_S(:);
		
		A        = sparse(  [id_S; id_S] ...
				  , [id_Sr(:); id_Sr(:)+1] ...
                                  , [p(:); 1-p(:)], nn,nr+1);
		A(:,end) = [];
%		Smu_  = A \ Shat(:);
	end

	% TODO interpolate for sd to
	%se = accumarray(flat(r),flat(f),[],@(x) std(x)/sqrt(length(x)));
end % periodogram_radial

