% Fri 18 Dec 19:19:50 +08 2020
function [Ni,si,vi] = oversampleNZ(N,sigma,v,d)
	nn = length(N);
	nz = length(sigma);
	s_ = inner2outer(cvec(sigma));
	s_ = min(max(s_,0),1);
	si = interp1(0:nz,s_,0:d:nz);
	si = mid(si);
	Ni = interp1(0:nn-1,N,0:d:nn-1);
	v = permute(v,[4,2,3,1]);
	vi = [];
	vi = interp1(-sigma,v,-si,'spline');
	%for cdx=1:size(z,1)
	%	vi(:,cdx,:,:) = interp1(-z(cdx,:),v(:,cdx,:,:),-zi(cdx,:),'spline'); %linear','extrap');
	%end
	vi = permute(vi,[4,2,3,1]);
	v  = permute(vi,[2,1,3,4]);
	v  = interp1(mid(N),v,mid(Ni),'spline');
	vi  = permute(v,[2,1,3,4]);
end

