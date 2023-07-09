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
%% analyze a 2D spatial pattern, estimate regularity and test for periodicity
%
function obj = analyze_grid(obj)
	obj.prepare_analysis();
	
	timer = tic();

	dx        = obj.stat.L_square ./ obj.stat.n_square;
	df_square = 1./obj.stat.L_square;

	% output
	S    = obj.S;
	R    = obj.R;
	stat = obj.stat;

	% statistics of the masked area
	area_msk = sum(obj.msk.b_square(:))*dx(1)*dx(2);
	% extend approximation by an ellipsis
	nmsk = sum(obj.msk.b_square(:));
	centroid.x = sum(obj.msk.b_square*cvec(obj.x))./nmsk;
	centroid.y = sum(rvec(obj.y)*obj.msk.b_square)./nmsk;
	sx2 = sum(obj.msk.b_square*(cvec(obj.x)-centroid.x).^2)./nmsk;
	sy2 = sum(obj.msk.b_square*(cvec(obj.y)-centroid.y).^2)./nmsk;
	sxy = sum(obj.msk.b_square*((cvec(obj.x)-centroid.x).*(cvec(obj.y)-centroid.y)))./nmsk;
	centroid.C = [sx2,sxy;
        	      sxy,sy2];

	% removal of spurious low frequency components
	% we assume that spurious low frequency components are predominantly isotropic
	[S.hp, whp, fhp, shp] = suppress_low_frequency_lobe(S.hat,obj.msk.b_square,obj.stat.L_square);
	%[S.hp, w] = suppress_low_frequency_components_1(obj,S.hat)

	% fraction of spectral energy retained after removing spurious low-frequency components by highpass filtering
	stat.p_S_hp = sum(whp.*S.hat,'all')./sum(S.hat,'all');

	% note that due to the high pass filtering, b_hp has positive and negative values
	b_hp = ifft2(sqrt(whp).*fft2(obj.b_square));

	% prescale image for thresholding
	% there is a matlab bug that double images need to be scaled
	b_scaled = (b_hp - min(b_hp(obj.msk.b_square)))/(max(b_hp(obj.msk.b_square))-min(b_hp(obj.msk.b_square)));
	stat.thresh_b      = graythresh(b_scaled(obj.msk.b_square));

	% threshold
	b_thresh      = b_scaled>stat.thresh_b;
	% fraction of ground coverged
	stat.coverage = sum(b_thresh(obj.msk.b_square(:)))/sum(obj.msk.b_square(:));

	% contrast between vegetated and unvegetated areas
	% computed for the original image as a quaility indicator
	stat.contrast = mean(obj.b_square(obj.msk.b_square & b_thresh)) - mean(obj.b_square(obj.msk.b_square & (~b_thresh)));

	obj.w.x = 1-normpdf(obj.f.x,0,shp)/normpdf(0,0,shp);	
	obj.w.y = 1-normpdf(obj.f.y,0,shp)/normpdf(0,0,shp);	
	obj.w.r = 1-normpdf(obj.f.r,0,shp)/normpdf(0,0,shp);

	% 2D spectral density estimate by periodogram smoothing
	dfr   = obj.f.rr(1,2);
	nf    = round(sqrt(obj.stat.q.fr.p50/dfr));
	S.bar = gaussfilt2(S.hp,nf);

	% smoothing window radius for frequency test
	nf_test  = sqrt(obj.opt.nf_test_scale*obj.stat.q.fr.p50/dfr);
	%nf_test  = round(0.25*obj.stat.q.fr.p50/dfr);
	nf_test  = max(nf_test,obj.opt.nf_test_min);

	% restrict test to region containing the upper 80% of spectral energy
	[Ssort,sdx] = sort(S.bar(:),'descend');
	iSsort = cumsum(Ssort);
	fdx = (iSsort <= obj.opt.p_fmsk*iSsort(end));
 	obj.msk.f = false(obj.stat.n_square);
	obj.msk.f(sdx(fdx)) = true;

	% exlude spurious low-frequency components from the test
	obj.msk.f = obj.msk.f & (obj.f.rr > fhp);

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
		[isperiodic, stati] = periodogram_test_periodicity_2d(...
					obj.b_square, obj.stat.L_square, nf_test, msk.b_, obj.msk.f_pos, obj.opt.n_mc);
		% note, for summing up the energy, the negative half-plane must not be excluded
		fdx                   = (stati.pn_all<=obj.opt.significance_level_a1) & obj.msk.f;
		p_periodic            = stati.pn;
		% as S is normalized to 1 over the half plane, this is identical to the fraction of spectral
		% energy contained in significant frequency components 
		stati.intS_hp_sig     = 0.5*sum(S.hp(fdx))*df_square(1)*df_square(2);
	else
		p_periodic = NaN;
		stati = struct();
	end
	%fr_periodic = stati.fr_max;
	catch e
		disp(e);
		p_periodic= NaN;
		stati = struct();
		%fr_periodic = NaN;
	end
	stat.stati = stati;

	% determine if pattern is isotropic (spotted, gapped, labyrinthic)
	% or anisotropic (banded)
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

	L_eff = effective_mask_size(obj.msk.b_square,obj.stat.L_square,-angle_deg);

	% for patterns with known direction, such as computer generated patterns, the direction angle_deg can be specified
	if (isfield('angle',obj.opt) && ~isempty(obj.opt.angle))
		angle_deg = obj.opt.angle;
	%else
		%angle = rad2deg(atan2(fc.yy.bar,fc.xx.bar));
	end

	L = obj.stat.L_square;
	df = 1./obj.stat.L_square;
	n  = obj.stat.n_square;
	for field = {'hp','hat','bar'}
		% normalize volume of density to 1
		S.(field{1})   = periodogram_normalize_2d(S.(field{1}),df_square);
		%S.(field{1})   = 2*S.(field{1})/(sum(sum(S.(field{1})))*df(1)*df(2));
		% autocorrelaton
		R.(field{1})   = 0.5*(n(1)*n(2))./(L(1)*L(2))*real(ifft2(S.(field{1})));


		% maximum of the 2d periodogram / density
		[Sc.(field{1}),cdx] = max(S.(field{1}),[],'all');
		%[stat.(field{1}).imax,stat.(field{1}).jmax] = ind2sub(n,cdx);
		[imax,jmax]    = ind2sub(n,cdx);
		% maximum frequency
		fc.rr.(field{1}) = obj.f.rr(cdx);
		fc.tt.(field{1}) = obj.f.tt(cdx);
		fc.xx.(field{1}) = obj.f.x(imax);
		fc.yy.(field{1}) = obj.f.y(jmax);

		% rotate periodogram/density
		S.rot.(field{1}) = fft_rotate(S.(field{1}),-angle_deg);
		R.rot.(field{1}) = fft_rotate(R.(field{1}),-angle_deg);

		% radial periodogram
		[S_,~] = periodogram_radial(S.(field{1}),L);
		S.radial.(field{1}) = S_.normalized;
	
		% angular periodogram
		[S.rot.angular.(field{1}),obj.f.angle] = periodogram_angular(S.rot.(field{1}),L,nf_s_); 


		% density perpendicular to bands
		S.rot.x.(field{1}) = mean(S.rot.(field{1}),2);
		% density parallel to bands
		S.rot.y.(field{1}) = mean(S.rot.(field{1}),1)';

		% normalize the area of density to 1 over the positive half-axis
		S.rot.x.(field{1}) = periodogram_normalize(S.rot.x.(field{1}),df_square(1));
		% normalize area under density to 1 over the entire axis
		S.rot.y.(field{1}) = periodogram_normalize(S.rot.y.(field{1}),df_square(2));

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
		[Sc.radial.(field{1}),id]  = max(S.radial.(field{1}));
		fc.radial.(field{1})       = obj.f.r(id);
		[Sc.x.(field{1}),id]       = max(cvec(S.rot.x.(field{1})).*(obj.f.x>=0));
		fc.x.(field{1})            = obj.f.x(id);
		[Sc.y.(field{1}),id]       = max(cvec(S.rot.y.(field{1})).*(obj.f.y>=0));
		% this should be the first field for anisotropic patterns, but we compute it anyway
		fc.y.(field{1})            = obj.f.y(id);
		[Sc.angular.(field{1}),id] = max(S.rot.angular.(field{1}));
		fc.angular.(field{1})      = obj.f.angle(id);

	end % for field

	obj.msk.rot.f = fft_rotate(obj.msk.f,-angle_deg);

	if (isisotropic)
		regularity = Sc.radial.hp .* fc.radial.hp;
		%[Sc.ref,ldx] = S.radial.hp);
		%lc       = 1./obj.f.r(ldx);
	else
		regularity = Sc.x.hp .* fc.x.hp;
		%[Sc,ldx] = max(Sx);
		%lc       = 1./abs(obj.f.x(ldx));
	end

	% store reasults
	stat.L_eff         = L_eff;
	stat.Sc	           = Sc;
	stat.angle_deg     = angle_deg;
	stat.area_msk      = area_msk;
	stat.fc            = fc;
	%stat.fr_periodic   = fr_periodic;
	stat.isisotropic   = isisotropic;
	stat.fhp           = fhp;
	stat.nf            = nf;
	stat.nf_test       = nf_test;
	stat.p_isotropic   = p_isotropic;
	stat.p_periodic    = p_periodic;
	stat.regularity    = regularity;
	stat.siz           = n;
	stat.centroid      = centroid;

	obj.stat = stat;
	obj.R = R;
	obj.S = S;

	obj.stat.analyzed = true;
	obj.stat = setfield_deep(obj.stat,'runtime.analysis',toc(timer));
end % analyze_grid

