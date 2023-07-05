% 2021-06-21 15:53:50.278752860 +0200
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
%%function [Sr,fri,se,count] = periodogram_radial(S2d,L)
%%
%% compute the radially averaged density
%%
%% input:
%% S2d : 2-dimensional density or periodogram
%% L =[Lx,Ly] : domain length
%%
%% output:
%%
%% S_r.mu         : radially averaged periodogram
%% S_r.normalized : normalized radially averaged periodogram
%% A              : matris operator s that Sr = (A*A')^-1 A'*S2d
%%
%% f_r : radial freqeuncies, at which radial periodogram is determined
%%       discretized in same interval as the 2d-density : f = 1/L
%%
%% Definitions:
%%        radial wavenumber, identical to circumferences of circles centred at origin with radial frequency fr
%%	       k_r = 2*pi*f_r
%%
%% 	 radially averaged periodogram:
%%        S_r(k_r)  = 1/k_r int_0^{k_r} S2d(k_r,s) d s
%%                  = 1/(2 pi) int_0^{2 pi} S2d(k_r,theta) d theta
%%                  ~ 1/(2 pi) sum^nt S2d(k_r,theta) * (2*pi/nt)
%%                  ~ 1/nt sum^nt S2d(k_r,theta)
%%                   
%%                nt ~ k_r/df = k_r*L
%% normalization:
%%       S_r.normalize  = S_r/int_0^inf S_r dfr
%%                      ~ S_r/(sum_0^nr S_r Delta fr)
%%
%% note : the radially averaged "periodogram", is actually a density estimate,
%%        for radial frequencies fr hat are not small
%%
%% when S is flattened into a vector, the isotropic part of the 2D density can be recovered with:
%% S_iso     = (A*S_radial)
%% S_radial  = A^-1 S_hat
%%
function [Sr, fri, count, A] = periodogram_radial(S2d,L)
	n = size(S2d);
	if (nargin()<2)
		L = n;
	end
	S     = struct();
	% input axes
	fx    = fourier_axis(L(1),n(1));
	fy    = fourier_axis(L(2),n(2));
	fr    = hypot(fx,fy');
	df    = 1./L;
	% note: this does not limit extend of radial density to circles which are entirely located inside of the rectangular domain
	frmax = hypot(max(fx),max(fy));

	% output axis
	dfi     = sqrt(df(1).*df(2));
%ceil(frmax+dfi))';
	fri_max = min(max(fx),max(fy));
	fri     = (0:dfi:fri_max);

	% (bi-linear) interpolation
	% each value is proportionally split between the next lower and next larger bin
	% bin index
	rat = fr/dfi;
	% integer part, index into 1D radial density
	id_Sr = floor(rat);
	% fractional part
	p   = 1-(rat-id_Sr);
	% account for 1-based indices in matlab, bin of fr=0 has index 1
	id_Sr = id_Sr+1;
	s = [length(fri)+1,1];

	fdx = id_Sr < s(1);

	% half-sum of 2d-bins per 1d-bin
if (0)
	sumS = 0.5*( accumarray(id_Sr(:),p(:).*S2d(:),s,@sum) ...
	           + accumarray(id_Sr(:)+1,(1-p(:)).*S2d(:), s, @sum) );
	% nt : number of 2d-bins per 1d-bin
	count = (   accumarray(id_Sr(:),p(:),s,@sum) ...
	          + accumarray(id_Sr(:)+1,(1-p(:)), s, @sum) );
else
	sumS = 0.5*( accumarray(id_Sr(fdx),p(fdx).*S2d(fdx),s,@sum) ...
	           + accumarray(id_Sr(fdx)+1,(1-p(fdx)).*S2d(fdx), s, @sum) );
	% nt : number of 2d-bins per 1d-bin
	count = (   accumarray(id_Sr(fdx),p(fdx),s,@sum) ...
	          + accumarray(id_Sr(fdx)+1,(1-p(fdx)), s, @sum) );
end

	% strip successor of last bin
	sumS  = sumS(1:end-1);
	count = count(1:end-1);

	% integral
	Sr.mu = sumS./count;
	Sr.mu(count==0) = 0;
	
	% normalize
	Sr.normalized = Sr.mu./(sum(mid(Sr.mu)).*dfi);

	% matrix
	if (nargout()>3)
		nr = length(fri);
		nn = numel(S2d);
		% index into 2D periodogram
		id_S     = (1:nn)';
		
		A        = sparse(  [id_S(fdx); id_S(fdx)] ...
				  , [id_Sr(fdx); id_Sr(fdx)+1] ...
                                  , [p(fdx); 1-p(fdx)], nn,nr+1);
		% strip successor of last bin
		A(:,end) = [];
		%Smu_  = A \ S2d(:);
	end

end % periodogram_radial

