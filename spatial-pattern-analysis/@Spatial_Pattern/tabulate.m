% 2021-07-06 15:17:53.166366045
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%% summarize properties of multiple patterns in a single struct
%% TODO function is deprecated, either remove or update
%
function tab = tabulate(obj)
	tab = struct();

	% tabulate properties
	tab.biomass               = arrayfun(@(x) x.stat.biomass,obj);
	tab.f.mean                = arrayfun(@(x) x.stat.f.mean, obj);
	tab.f.mean0               = arrayfun(@(x) x.stat.f.mean0, obj);
	tab.f.std                 = arrayfun(@(x) x.stat.f.std, obj);
	tab.f.skew                = arrayfun(@(x) x.stat.f.skew, obj);
	tab.f.kurt                = arrayfun(@(x) x.stat.f.kurt, obj);
	f_C = {'bandpass','brownian','lognormal','lorentzian'};
	for idx=1:length(f_C)
	    tab.good.(f_C{idx})       = arrayfun(@(x) x.stat.good.(f_C{idx}), obj);
	    tab.par1.(f_C{idx})       = arrayfun(@(x) x.stat.par1.(f_C{idx}), obj);
	    tab.par2.(f_C{idx})       = arrayfun(@(x) x.stat.par2.(f_C{idx}), obj);
	    tab.r2.(f_C{idx})         = arrayfun(@(x) x.stat.r2.(f_C{idx}), obj);
	    tab.Sc.(f_C{idx})         = arrayfun(@(x) x.stat.Sc.(f_C{idx}), obj);
	    tab.wavelength.(f_C{idx}) = arrayfun(@(x) x.stat.wavelength.(f_C{idx}), obj);
	end
	tab.good.bartlett         = arrayfun(@(x) x.stat.good.bartlett, obj);
	f_C = {'gauss','bartlett','filt'};
	for idx=1:length(f_C)
		tab.Sc.(f_C{idx})               = arrayfun(@(x) x.stat.Sc.(f_C{idx}), obj);
		tab.wavelength.(f_C{idx})       = arrayfun(@(x) x.stat.wavelength.(f_C{idx}),obj);
%		tab.Sc.filt               = arrayfun(@(x) x.stat.Sc.bartlett, obj);
%		tab.wavelength.filt       = arrayfun(@(x) x.stat.wavelength.filt,obj);
%		tab.Sc.bartlett           = arrayfun(@(x) x.stat.Sc.bartlett, obj);
%		tab.wavelength.bartlett   = arrayfun(@(x) x.stat.wavelength.bartlett,obj);
	end

	tab.p_periodic            = arrayfun(@(x) x.stat.p_periodic, obj);
	tab.periodic              = arrayfun(@(x) x.stat.periodic, obj);

	try % only defined for modelled patterns
	tab.celerity              = arrayfun(@(x) x.stat.celerity,obj);
	catch e
	end
	try
	tab.corr		  = arrayfun(@(x) x.stat.corr,obj);
	catch e
	end
	tab.length_m              = arrayfun(@(x) x.L,obj);
	tab.nx                    = arrayfun(@(x) x.n,obj);
	tab.mseg                  = arrayfun(@(x) x.stat.m, obj);

end

