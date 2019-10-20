% Fri  5 May 13:31:18 CEST 2017
%
%% class for short time fourier transform
%%
%% Note: the interval Ti should be set to at leat 2*max(T), as otherwise coefficients
%%       tend to oscillate in the presence of noise
%% Note: for convenience, the independent variable is labeled as time (t),
%%       but the independent variable is arbitrary, so it works likewise in space
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
classdef STFT < handle
	properties
		% frequencies to consider
		T
		% winfun

		% mesh quantities		

		% start time
		t0
		% end time
		tend

		% interval length
		Ti

		% order of interpolation polynomial
		order = 3;
	
		% fourier coefficients
		coeff

		% discretisation matrices
		A1
		A2
		% maximum iteration of minres
		maxit  = 1000;
		lambda = 0;
	end % properties
	methods
		function obj = STFT(varargin)
			for idx=1:2:length(varargin)-1
				field = varargin{idx};
				val   = varargin{idx+1};
				obj.(field) = val;	
			end % for
			% quick fix
			if (3 == obj.order)
				obj.t0 = obj.t0-obj.Ti;
	
			end
		end % constructor
	
		% pseudo member variables (just in time computation)
		function [ic obj] = icoeff(obj,t)
			if (nargin()>1)
				[void coeff] = obj.itransform(t);
			else
				coeff = obj.coeff;
			end
			ic            = zeros(1+length(obj.T),size(coeff,2),(size(coeff,3)+1)/2);
			ic(1,:,:)     = coeff(1,:,:);
			% TODO swap im and re?
			ic(2:end,:,:) = coeff(2:2:end-1,:,:)+1i*coeff(3:2:end,:,:);
		end

		function [a obj] = amplitude(obj,t)
			if (nargin()>1)
				[void coeff] = obj.itransform(t);
			else
				coeff = obj.coeff;
			end
			a            = zeros(1+length(obj.T),size(coeff,2),size(coeff,3));
			a(1,:,:)     = abs(coeff(1,:,:));
			a(2:end,:,:) = hypot(coeff(2:2:end-1,:,:),coeff(3:2:end,:,:));
		end
		function [p obj] = phase(obj,t)
			if (nargin()>1)
				[void coeff] = obj.itransform(t);
			else
				coeff = obj.coeff;
			end
			p            = zeros(1+length(obj.T),size(coeff,2),size(coeff,3));
			p(1,:,:)     = angle(coeff(1,:,:)); % sign
			p(2:end,:,:) = atan2(coeff(3:2:end,:,:),coeff(2:2:end-1,:,:));
		end
		function obj = set_phase(obj,p,a)
			if (nargin()<3)
				a = obj.amplitude();
			end
			c = cos(p);
			s = sin(p);
			% do not use end, as coeff can be empty
			ni = obj.ni;
			obj.coeff(1,:,:)        = sign(c(1,:)).*a(1,:,:);
			obj.coeff(2:2:ni-1,:,:) = c(2:end,:).*a(2:end,:,:);
			obj.coeff(3:2:ni,:,:)   = s(2:end,:).*a(2:end,:,:);
			%obj.coeff = c.*a;
			%obj.coeff = s.*a;
		end
		function obj = set_amplitude(obj,a,p)
			if (nargin()<3)
				p = obj.phase;
			end
			c = cos(p);
			s = sin(p);
			% do not use end, as coeff can be empty
			ni = obj.ni;
			obj.coeff = c.*a;
			obj.coeff = s.*a;
			%obj.coeff(1,:,:)        = -a(1,:,:);
			obj.coeff(1,:,:)        = sign(c(1,:)).*a(1,:,:);
			obj.coeff(2:2:ni-1,:,:) = c(2:end,:).*a(2:end,:,:);
			obj.coeff(3:2:ni,:,:)   = s(2:end,:).*a(2:end,:,:);
		end
		function [ni, obj] = ni(obj)
			ni = 1+2*length(obj.T);
		end
		function [nti, obj] = nti(obj)
			nti = floor((obj.tend-obj.t0)/obj.Ti)+2;
		end
		% time of supporting coefficients
		function [ti  obj] = ti(obj)
			ti = obj.t0 + obj.Ti*(0:obj.nti-1)';
		end
	end % methods
	methods (Static)
		[tc c amp phase serr] = stft(val,dt,Ti,T,winstr);
	end
end % classdef STFT

