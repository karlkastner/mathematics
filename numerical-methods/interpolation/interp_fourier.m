% Di 26. Jan 12:32:46 CET 2016
% Karl Kastner, Berlin
%
%% interpolation by the fourier method
function [tF serr] = inter_fourier(sS,sNrel,sF,tS,tN)
	nt = size(tS);
	ns = size(sS);
	serr = NaN;
	tF = NaN(size(tS));
	if (nt(1) < 2 || nt(2) < 2)
		return;
	end

	% number of eigenfunctions in series expansion
%	k = round(sqrt(numel(tS)));
	k = ceil(numel(tS)/100);

	% TODO estimate w from data
	w = nanmedian(abs(tN(:,1)-tN(:,end)));
	L = [abs(tS(end,end)-tS(1,1)),w];

	% TODO, the eigenfunctions are just cosines, generate automatically
	% set up Laplacian operator
	A = laplacian_fdm(nt,L,'neumann');

	% determine the lowest frequent eigenfunctions
	% of the Laplacian at the grid points
	opts.issym = true;
	%[V E] = eigs(A,k,'SM',opts);
	[V E] = eigs(A,k,sqrt(eps),opts);

	% TODO, as the eigenfunctions are just cosines, they can be expanded and do not need to
	% be interpolated

	% interpolate eigenfunctions to sample locations
	tNrel = (1:nt(2))/(nt(2)+1)-0.5;
	Vi = zeros(numel(sS),k);
	for idx=1:k
%		V_ = reshape(V(:,idx),nt(1),nt(2));
		V_ = reshape(V(:,idx),nt(2),nt(1));
%		clf
%		imagesc(V_)
%		pause
		Vi(:,idx) = interp2(tS(:,1),tNrel,V_,flat(sS),flat(sNrel));
	end
	fdx = isfinite(sum(Vi,2));
	if (sum(fdx) < k+1)
		return;
	end

	% determine coefficients by solving Vi c = F
	sF = flat(sF);
	%[c Rcond] = linsolve(Vi(fdx,:),sF(fdx));
	[c] = robustfit(Vi(fdx,:),sF(fdx),[],[],false);
	
	% estimate the error from the residual
	res  = Vi*c - sF;
	ns_ = sum(fdx);
	serr = sqrt(res(fdx)'*res(fdx)/(ns_-k));

	% expand function at grid points as series (weighed sum) of eigenfunctions
	tF = V*c;
	tF = reshape(tF,nt(2),nt(1))';
if (0)
	clf
	%imagesc(tS(:,1),tN(1,:),tF)
	%surface(tS(:,1),tN(1,:),[],tF)
	%imagesc(tN(1,:),tS(:,1),tF)
	surface(tS(:,1)*ones(1,nt(2)),ones(nt(1),1)*tN(1,:),tF,'edgecolor','none')
	ylim([-0.5 0.5])
	hold on
	plot(sS,sNrel,'r.')
	colorbar
	pause
end
end

