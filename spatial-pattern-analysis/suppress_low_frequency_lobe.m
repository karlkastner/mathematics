% 2022-12-02 00:25:33.449376625 +0100
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
function [Shat, whp, fhp, shp, nf] = suppress_low_frequency_lobe(Shat,msk,L)
	siz = size(Shat);
	[fx_i,fy_i,frr_i] = fourier_axis_2d([1,1],siz);
	[fx,fy,frr] = fourier_axis_2d(L,siz);

	% radial density, including low-frequency components
	[Sr, fr] = periodogram_radial(Shat,L);
	Sr0 = Sr.normalized;
	% radial distribution (cdf)
	iS = cumsum(Sr.normalized)*(fr(2)-fr(1));
	iS = iS/iS(end);

	% smoothing window size ensuring convergence	
	nf = round(sqrt(find(iS>0.5,1,'first')));

	% provisional estimate of the density including the spurious low-frequency
	% components by smoothing the periodogram
	Shat_smooth = Shat;
	%Shat_smooth(1,1) = 0.5*(Shat(1,2)+Shat(2,1));
	L_eff = effective_mask_size(msk,L,0);

	l_ = round(sqrt(L(1)*L(2))./L_eff.r);
	%l_ = ceil(sqrt(numel(msk)./sum(msk(:))))+1;
	% the zero at the origin has to be filled for the min-max algorithm to work
	% masking can smooth the zero across several pixels, so several pixels have to be filled
	Shat_smooth(frr_i <= l_) = max(Shat_smooth(frr_i<=l_+1));
	%Shat_smooth(frr_i <= l_) = max(Shat_smooth(frr_i>l_ & frr_i<=l_+1));
%	else
		% proxy for mean (4 point averge, reduced to 2 points due to symmetry)
		
%	end
	Sr_ = periodogram_radial(Shat_smooth);

	% TODO gaussian smoothing window
	%Shat_smooth = circfilt2(Shat_smooth,nf);
	Shat_smooth = gaussfilt2(Shat_smooth,nf);
	%Shat_smooth = ifftshift(trifilt2(fftshift(Shat_smooth),nf));

	% re-zero mean
	Shat_smooth(1,1) = 0;

	% radial periodogram of smoothed periodogram
	Sr_smooth = periodogram_radial(Shat_smooth,L);

	% minmax of spectral energy,
	% the frequency of the minimum is used for separating the densities
	[mindx,maxdx] = minmax(Sr_smooth.normalized(2:end));
	mindx = mindx+1;	
	fhp = fr(mindx);

	% remove spurious low-frequency components by multiplication
	% with Fourier window
	if(0)
		dfrr = hypot(fx(2)-fx(1),fy(2)-fy(1));
		shp = (fhp)/sqrt(2*log(2));
		wf = 1-normpdf(frr,0,shp)/normpdf(0,0,shp);
	else
		%dfrr = hypot(fx(2)-fx(1),fy(2)-fy(1));
		shp  = (fhp)/sqrt(2*log(2));
		whp  = 1-normpdf(frr,0,shp)/normpdf(0,0,shp);	
	end

	%Shat(frr_i<=hp(idx)) = 0;
	%s = (hp(idx))/sqrt(2*log(2))
	%w = 1-normpdf(frr_i,0,s)/norm;
	Shat = whp.*Shat;

end

