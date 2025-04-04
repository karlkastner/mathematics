% Wed 18 May 13:50:47 CEST 2022
% Karl Kästner, Berlin
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%% analyze a spatial pattern in two-dimensions, estimate regularity and
%% test for periodicity
%
function obj = analyze_grid(obj)
	obj.prepare_analysis();
	
	timer = tic();

	% spectral resolution of pattern extended to square dimensions
	df_square = 1./obj.stat.L_square;

	% for quicker access, R,S,stat are fetched here and later written back
	S    = obj.S;
	T    = obj.T;
	R    = obj.R;
	stat = obj.stat;

	% statistics of the masked area
	area_msk = obj.area_msk();
	centroid = obj.centroid();	

	% removal of spurious low frequency components
	% assumes that spurious low frequency components are predominantly isotropic
	if (obj.opt.suppress_low_frequency_components)
	[S.hp, whp, f_hp, sd_hp, nf_lf, dSminmaxrel] = suppress_low_frequency_lobe(S.hat,obj.msk.b_square,obj.stat.L_square);
	else
		S.hp = S.hat;
		whp = 1;
		f_hp = 0;
		dSminmaxrel = NaN;
	end
	%[S.hp, w] = suppress_low_frequency_components_1(obj,S.hat)

	% fraction of spectral energy retained after removing spurious low-frequency components by highpass filtering
	stat.p_S_hp = sum(whp.*S.hat,'all')./sum(S.hat,'all');

	% note that due to the high pass filtering, b_hp has positive and negative values
	b_hp = ifft2(sqrt(whp).*fft2(obj.b_.square));

	% prescale image for thresholding
	% there is a matlab bug that double images need to be scaled
	b_scaled = (   (b_hp - min(b_hp(obj.msk.b_square))) ...
		     / (max(b_hp(obj.msk.b_square)) ...
		     -  min(b_hp(obj.msk.b_square))));
	stat.thresh_b      = graythresh(b_scaled(obj.msk.b_square));

	% threshold
	obj.b_.thresh      = b_scaled>stat.thresh_b;

	% fraction of ground covered by vegetation
	stat.coverage = (   sum(obj.b_.thresh(obj.msk.b_square),'all') ...
		          / sum(obj.msk.b_square,'all') );

	% contrast between vegetated and unvegetated areas
	% computed for the original image as a quaility indicator
	stat.contrast = (  mean(obj.b_.square(obj.msk.b_square & obj.b_.thresh)) ...
			 - mean(obj.b_.square(obj.msk.b_square & (~obj.b_.thresh))));

	% weights for density fits, exclude spurious low-frequency components
	if (obj.opt.suppress_low_frequency_components & obj.opt.weight)
		obj.w.x = 1-normpdf(obj.f.x,0,sd_hp)/normpdf(0,0,sd_hp);	
		obj.w.y = 1-normpdf(obj.f.y,0,sd_hp)/normpdf(0,0,sd_hp);	
		obj.w.r = 1-normpdf(obj.f.r,0,sd_hp)/normpdf(0,0,sd_hp);
	else
		obj.w.x = ones(size(obj.f.x));
		obj.w.x(1) = 0;
		obj.w.y = ones(size(obj.f.y));
		obj.w.r = ones(size(obj.f.r));
		obj.w.r(1) = 0;
	end
	obj.w.xp = obj.w.x.*(obj.f.x>0);

	% 2D spectral density estimate by periodogram smoothing
	dfr   = obj.f.rr(1,2);
	nf    = round(sqrt(obj.stat.q.fr.p50/dfr));
	S.bar = gaussfilt2(S.hp,nf);

	% smoothing window radius for frequency test
	nf_test  = sqrt(obj.opt.nf_test_scale*obj.stat.q.fr.p50/dfr);
	nf_test  = max(nf_test,obj.opt.nf_test_min);

	% restrict test to region containing the upper 80% of spectral energy
	[Ssort,sdx] = sort(S.bar(:),'descend');
	iSsort = cumsum(Ssort);
	fdx = (iSsort <= obj.opt.p_fmsk*iSsort(end));
 	obj.msk.f = false(obj.stat.n_square);
	obj.msk.f(sdx(fdx)) = true;

	% exlude spurious low-frequency components from the test
	obj.msk.f = obj.msk.f & (obj.f.rr > f_hp);

	% by symmetry, the negative-half plane can be excluded
	obj.msk.f_pos = obj.msk.f & cvec(obj.f.x) >= 0;

	% periodicity test
	try
	if (any((obj.msk.b_square~=obj.msk.b_square(1,1)),'all'))
		msk.b_ = obj.msk.b_square;
	else
		msk.b_ = [];
	end
	if (obj.opt.test_for_periodicity)
		[isperiodic, p_periodic, stati] = periodogram_test_periodicity_2d(...
					obj.b_.square, obj.stat.L_square, nf_test, msk.b_, obj.msk.f_pos, obj.opt.n_mc);
		% note, for summing up the spectral energy, the negative half-plane must not be excluded
		fdx                   = (stati.pn_all<=obj.opt.significance_level_a1) & obj.msk.f;
		% p_periodic            = stati.pn;
		% as S is normalized to 1 over the half plane, this is identical to the fraction of spectral
		% energy contained in significant frequency components 
		stati.intS_hp_sig     = 0.5*sum(S.hp(fdx))*df_square(1)*df_square(2);
	else
		p_periodic = NaN;
		stati = struct();
	end
	catch e
		disp(e);
		disp(e.stack);
		p_periodic = NaN;
		stati = struct('e',e);
	end
	stat.stati = stati;

	% determine if pattern is isotropic (spotted, gapped, labyrinthic)
	% or anisotropic (striped)
	% note: presmoothed with Sbar works better than with Shat
	%nt   = pi*n(1);
	%nc   = 2*pi*f_50/dfr;
	nc   = sqrt(obj.stat.q.fr.p50/dfr);
	nf_s = ceil(obj.stat.n_square(1)/nc);
	% more sharper later
	nf_s_  = ceil(1/4*obj.stat.n_square(1)/nc);
	mode   = 'angular'; 
	nf_    = round(2*sqrt((obj.stat.q.fr.p50/dfr)));
	Sbar_  = gaussfilt2(S.hp,nf_);
	[isisotropic,stati] = separate_isotropic_from_anisotropic_density(Sbar_,obj.msk.f,obj.stat.L_square,mode,nf_s);
	angle_deg = stati.angle_deg;
	p_isotropic = stati.p_iso;

	% effective spatial extent of masked area
	L_eff = effective_mask_size(obj.msk.b_square,obj.stat.L_square,-angle_deg);

	% smoothing window for consitently estimating the spectral density 
	% of patterns with spatial extent of high aspect ratio
	% target spectral resolution
	Lt = sqrt(L_eff.x*L_eff.y);
	n = size(S.hat);
	if (obj.stat.L_square(1) > Lt)
		% when the pattern is shorter along the x-axis, then the
		% y-component is averaged over fewer samples and has to be smoothed
		%m   = sqrt(L_eff.x/L_eff.y);
		nwx  = (obj.stat.L_square(1)/Lt);
		swx  = gausswin_dof2std(nwx,obj.f.x(2)-obj.f.x(1));
		% gaussian window equivalent to m-degrees of freedom
		winx = normpdf(obj.f.x,0,swx);
	else
		winx = zeros(n(1),1);
		winx(1) = 1;
	end
 	if (obj.stat.L_square(2) > Lt)
		% when the pattern is shorter along the y-axis, then the
		% x-component is averaged over fewer samples and has to be smoothed
		%m   = sqrt(L_eff.y./L_eff.x);
		nwy   = (obj.stat.L_square(2)./Lt);
		% gaussian window equivalent to m-degrees of freedom
		swy  = gausswin_dof2std(nwy,obj.f.y(2)-obj.f.y(1));
		winy = normpdf(obj.f.y,0,swy);
	else
		winy = zeros(1,n(2));
		winy(1) = 1;
	end
	Ws = cvec(winx)*rvec(winy);
	% the tails of the gaussian might be truncated, thus normalization
	% is necessary to ensure unit volume
	Ws = Ws/sum(Ws,'all');

	% for patterns with known direction, such as computer generated patterns, the direction angle_deg can be specified
	if (isfield('angle_deg',obj.opt) && ~isempty(obj.opt.angle_deg))
		angle_deg = obj.opt.angle;
	end

	L = obj.stat.L_square;
	df = 1./obj.stat.L_square;
	n  = obj.stat.n_square;
	% prepare output
	fc = struct();
	Sc = struct();
	regularity = struct();
	for field = {'hp','hat','bar'}
		% normalize volume of density to 1
		S.(field{1})   = periodogram_normalize_2d(S.(field{1}),df_square);
		%S.(field{1})   = 2*S.(field{1})/(sum(sum(S.(field{1})))*df(1)*df(2));
		% autocorrelaton
		R.(field{1})   = (n(1)*n(2))./(L(1)*L(2))*real(ifft2(S.(field{1})));

		% rotate periodogram/density
		S.rot.(field{1}) = fft_rotate(S.(field{1}),-angle_deg);
		R.rot.(field{1}) = fft_rotate(R.(field{1}),-angle_deg);

		% maximum of the 2d periodogram / density
		[Sc.(field{1}),cdx] = max(S.(field{1}),[],'all');
		%[stat.(field{1}).imax,stat.(field{1}).jmax] = ind2sub(n,cdx);
		[imax,jmax]    = ind2sub(n,cdx);
		% maximum frequency
		fc.rr.(field{1}) = obj.f.rr(cdx);
		fc.tt.(field{1}) = obj.f.tt(cdx);
		fc.xx.(field{1}) = obj.f.x(imax);
		fc.yy.(field{1}) = obj.f.y(jmax);

	end

	% convolve
	% TODO this can be optimized by two one-dimensional filtering operations
	% along each direction instead of one two-dimensional
	% TODO it would be better to fit an ellipse to the spatial extent and
	%      then to filter along the direction of the axes
	R.rot.con = fft2(Ws).*fft2(S.rot.hp);
	S.rot.con = ifft2(R.rot.con);
	R.rot.con = real(R.rot.con);
	R.rot.con = R.rot.con/R.rot.con(1);
	% rotate back
	S.con = fft_rotate(S.rot.con,+angle_deg);
	R.con = fft_rotate(R.rot.con,+angle_deg);

	if (~isempty(obj.source))
		R.e.rot.con  = fft2(Ws).*fft2(S.e.hat);
		S.e.rot.con  = ifft2(R.e.rot.con);
		R.be.rot.con = fft2(Ws).*fft2(S.be.hat);
		S.be.rot.con = ifft2(R.be.rot.con);
		obj.T.rot.con = S.be.rot.con/S.e.rot.con;
	end

	for field = {'hp','hat','bar','con'}
		% radial periodogram
		[S_,~] = periodogram_radial(S.(field{1}),L);
		S.radial.(field{1}) = S_.normalized;
		% angular periodogram
		[S.angular.(field{1}),obj.f.angle] = periodogram_angular(S.(field{1}),L,nf_s_);
	
		% angular periodogram, rotated
		[S.rot.angular.(field{1}),obj.f.angle] = periodogram_angular(S.rot.(field{1}),L,nf_s_);

		% over half-circle
		S.rot.angular_p.(field{1}) = ( 2*cvec(S.rot.angular.(field{1})) ...
			                      .* (  cvec(obj.f.angle) >= -pi/2 ...
			                          & cvec(obj.f.angle) <   pi/2 ...
                                                 ) ...
                                             );


		% density perpendicular to bands
		S.rot.x.(field{1})  = sum(S.rot.(field{1}),2)*df(2);

		% over positive half-axis
		S.rot.xp.(field{1}) = ( 2*S.rot.x.(field{1}) .* (cvec(obj.f.x) >= 0) );

		% density parallel to bands
		S.rot.y.(field{1}) = (sum(S.rot.(field{1}),1)')*df(1);

		% normalize the area of density to 1 over the positive half-axis
		%fullaxsis = false;
		%S.rot.x.(field{1}) = periodogram_normalize(S.rot.x.(field{1}),df_square(1),false);

		% normalize area under density to 1 over the entire axis
		%fullaxis = true;
		%S.rot.y.(field{1}) = periodogram_normalize(S.rot.y.(field{1}),df_square(2),true);

		% autocorrelation in direction perpendictular to bands
		R.rot.x.(field{1}) = mean(R.rot.(field{1}),2);
		R.rot.x.(field{1}) = R.rot.x.(field{1})/R.rot.x.(field{1})(1);
		%R.rot.x.(field{1}) = 0.5*n(1)/L(1)*real(ifft(S.rot.x.(field{1})));

		% autocorrelation in direction parallel to bands
		%R.rot.y.(field{1}) = 0.5*n(2)/L(1)*real(ifft(S.rot.y.(field{1})));
		R.rot.y.(field{1}) = mean(R.rot.(field{1}),1)';
		R.rot.y.(field{1}) = R.rot.y.(field{1})/R.rot.y.(field{1})(1);

		% radial autocorrelation
		[R.radial.(field{1}),obj.r] = autocorr_radial(R.(field{1}),L);

		% density maxima
		[Sc.radial.([field{1},'_']),id]  = max(S.radial.(field{1}));
		 fc.radial.([field{1},'_'])       = obj.f.r(id);
		[Sc.radial.(field{1}),fc.radial.(field{1})] = extreme3(obj.f.r,S.radial.(field{1}),id);
		% this is symmetric and we search here only in the positive half
		[Sc.x.([field{1},'_']),id]          = max(cvec(S.rot.x.(field{1}).*(obj.f.x>=0)));
		fc.x.([field{1},'_'])            = obj.f.x(id);
		[Sc.x.(field{1}),fc.x.(field{1})] = extreme3(obj.f.x,S.rot.x.(field{1}),id);
		[Sc.xp.([field{1},'_']),id]      = max(cvec(S.rot.xp.(field{1})));
		fc.xp.([field{1},'_'])            = obj.f.x(id);
		[Sc.xp.(field{1}),fc.xp.(field{1})] = extreme3(obj.f.x,S.rot.xp.(field{1}),id);
		%[Sc.xp.(field{1}),id]      = max(cvec(S.rot.xp.(field{1})));
		[Sc.y.(field{1}),id]       = max(cvec(S.rot.y.(field{1})).*(obj.f.y>=0));
		% this should be 0 for anisotropic patterns, but we compute it anyway
		fc.y.(field{1})            = obj.f.y(id);
		[Sc.angular.(field{1})]    = max(S.rot.angular.(field{1}));
		[Sc.angular_p.(field{1}),id] = max(S.rot.angular_p.(field{1}));
		fc.angular.(field{1})      = obj.f.angle(id);

		if (isisotropic)
			regularity.(field{1}) = Sc.radial.(field{1}) .* fc.radial.(field{1});
		else
			regularity.(field{1}) = Sc.xp.(field{1}) .* fc.x.(field{1});
		end
	end % for field

	if (~isempty(obj.source))
		S.e.rot.x.con  = sum(S.e.rot.con,2)*df(2);
		S.be.rot.x.con = sum(S.be.rot.con,2)*df(2);
		S.e.rot.y.con  = (sum(S.e.rot.con,1)')*df(1);
		S.be.rot.y.con = (sum(S.be.rot.con,1)')*df(1);
		T.rot.x = S.be.rot.x.con ./ S.e.rot.x.con;
		T.rot.y = S.be.rot.y.con ./ S.e.rot.y.con;
		% average radial coherence
		rSr = cvec(obj.f.r).*cvec(S.radial.hat);
		df = obj.f.r(2)-obj.f.r(1);
		rSr = rSr/(sum(rSr)*df);
		stat.coherence.radial = rSr'*obj.S.coherence.radial*df;
	end

	obj.msk.rot.f = fft_rotate(obj.msk.f,-angle_deg);

	% store reasults
	stat.L_eff         = L_eff;
	stat.Sc	           = Sc;
	stat.angle_deg     = angle_deg;
	stat.area_msk      = area_msk;
	stat.fc            = fc;
	stat.isisotropic   = isisotropic;
	stat.f_hp           = f_hp;
	stat.nf            = nf;
	stat.nf_test       = nf_test;
	stat.p_isotropic   = p_isotropic;
	stat.p_periodic    = p_periodic;
	stat.regularity    = regularity;
	stat.siz           = n;
	stat.centroid      = centroid;
	stat.dSminmaxrel   = dSminmaxrel;
	stat.Ws = Ws;

	obj.stat = stat;
	obj.R = R;
	obj.S = S;
	obj.T = T;

	obj.stat.analyzed = true;
	obj.stat = setfield_deep(obj.stat,'runtime.analysis',toc(timer));

end % analyze_grid

