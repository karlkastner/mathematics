% Mon  5 Dec 14:28:58 CET 2022
function [S,f,piso,slope] = separate_isotropic_from_anisotropic_density(Shat,L,nf)
	n = size(Shat);
	if (n(1)~=n(2) || L(2)~=L(1))
		error('periodogram must be square as it will distort by rotation')
	end

	S = struct();
	f = struct();
        [f.x,f.y,f.r] = fourier_axis_2d(L,size(Shat));                             
        fxx = repmat(cvec(f.x),1,length(f.y));                                    
        fyy = repmat(rvec(f.y),length(f.x),1); 


	mode = 'max';

	switch (mode)
	case {'ls'}
	% least squares
	A = [ones(numel(Shat),1),fxx(:)];
	w = Shat(:);
	W = diag(sparse(Shat(:)));
	c = (A'*W*A) \ (A'*W*fyy(:));
	slope = c(2);
	if (abs(c(2))>1)
		% swap to improve numerical accuracy
		A = [ones(numel(Shat),1),fyy(:)];
		%Shat_ = rot90(Shat);
		W = diag(sparse(Shat(:)));
		c = (A'*W*A) \ (A'*W*fxx(:));	
		slope = 1./c(2);
	end
	case {'po'}
		% this somethins yields incorrectly 0
	        slope = least_squares_perpendicular_offset(fxx(:),fyy(:),Shat(:));     
	case {'max'} 
		Shat_ = ifftshift(trifilt2(fftshift(Shat),nf));
		[Sc,mdx] = max(Shat_(:));
		if (fxx(mdx)==0)
			slope = inf;
		else
		slope = fyy(mdx)./fxx(mdx);
		end
	end

        a = atand(slope); 
        S.hat = ifftshift(imrotate(fftshift(Shat),-a,'crop','bilinear')); 

	% mask extra-diagonal pattern
	%x = zeros(100,100);
	%[fx,fy,fr] = fourier_axis_2d([1,1],[100,100]);
	S.msk = (abs(f.x)>abs(f.y'));
	%Shat_msk      = S.hat;
	%Shat_msk(~S.msk) = 0;
	S.hat_iso_msk = S.hat.*(~S.msk);
	% estimate radial density from excluded part
	[S.r,f.r,~,A] = periodogram_radial(S.hat.*(~S.msk),L);
	S.iso = reshape(A*S.r.normalized,n);
	% scale energy
	S.iso = (2*sum(flat(S.hat.*(~S.msk))))/sum(S.iso(:))*S.iso;
	S.hat_aniso = S.hat - S.iso;
	S.y = mean(S.hat_aniso)';
	S.x = mean(S.hat_aniso,2);
	S.x = S.x/sum(S.x(f.x>=0)*(f.x(2)-f.x(1)));
	S.y = S.y/sum(S.y(f.y>=0)*(f.y(2)-f.y(1)));
	S.aniso = S.x*S.y';
	% factor for anisotropic part
	% alternatively, we can raise S to the power
	if (0)
	fa = 2.5; % for 2, aniso are missclassified as iso, for 3, iso are misscl. as ansio
	piso = (sum(S.iso(:)) - fa*sum(S.hat_aniso(:)))./sum(Shat(:));
	else
		p = 1.5;
		e_iso   = 2*sum(flat(S.hat_iso_msk.^p));
		e_aniso =   sum(flat(max(S.hat_aniso,0).^p));
		piso = (e_iso - e_aniso)./(e_iso + e_aniso);
	end
end

