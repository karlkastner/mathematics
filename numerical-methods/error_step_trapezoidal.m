% Tue 14 Nov 13:01:46 CET 2023
% Karl Kastner, Berlin
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
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
% leading term of the local error
%% truncation error
%% note that the global error is larger!!!
function [e,dt0,out] = error_step_trapezoidal(t,z,erel0,final)
	if (nargin()<4 || ~final)
		% compute error for entire time series
		T  = t(end)-t(1);
		nt = length(t);
		dt = T/(nt-1);
		% TODO variable step
		D3 = derivative_matrix_k_1d(nt,T,3);
		d3z_dt3 = D3*double(z);
	else
		% expects last value in tt(1) and z(:,1)
		A        = vander_1d(t-t(1),3);
		v        = vanderd_1d(0,3,3);
		%d3z_dt3    = v*(z*inv(A)')';
		d3z_dt3    = (z*(A'\v'));
		dt = t(1)-t(2);
		%e        = -0.5*d2z_dt2*dt*dt;
		%e        = e_div_dt2*dt^2;
	end
	e    = (-dt*dt*dt/12)*d3z_dt3;
	rmsz = rms(z(:,1));
	% einf = max(abs(e),[],2);
	erms = rms(e);
	erel = erms./rmsz;
	if (nargin()>2)
		dt0 = dt*cbrt(erel0./abs(erel));
	else
		dt0 = [];
	end
	if (nargout()>2)
		%out.zrms = zrms;
		out.rmsz = rmsz;
		out.erms = erms;
		out.erel = erel;
		out.dt0  = dt0;
	end
end

