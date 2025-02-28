% 2022-09-26 14:32:10.449630621 +0200
% Karl Kästner, Berlin
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
function shp = export_shp(obj,shpname)

	mydeal              = @(x) deal(x{:});
	obj.runtime(end:obj.n,:) = NaN;

	% create shapefile
	n                 = size(obj.centroid,1);
	shp               = repmat(struct(),n,1);
	[shp.X]           = mydeal(num2cell(obj.centroid(:,1)));
	[shp.Y]           = mydeal(num2cell(obj.centroid(:,2)));
	[shp.L_eff_r]     = mydeal(num2cell(obj.L_eff_r));
	[shp.L_eff_x]     = mydeal(num2cell(obj.L_eff_x));
	[shp.area_m2]     = mydeal(num2cell(obj.area_msk));
	[shp.contrast]    = mydeal(num2cell(obj.contrast));
	[shp.coverage]    = mydeal(num2cell(obj.coverage));
	[shp.filename]    = mydeal(obj.filename_C);
	[shp.intS_hp_sig] = mydeal(num2cell(obj.intS_hp_sig));
	[shp.isisotropic] = mydeal(num2cell(obj.isisotropic));
	[shp.p_periodic]  = mydeal(num2cell(obj.p_periodic));
	[shp.p_S_hp]      = mydeal(num2cell(obj.p_S_hp));
	[shp.region_id]   = mydeal(num2cell(obj.region_id));
	[shp.region_str]  = mydeal(obj.region_C(obj.region_id));
	[shp.regulari_r]  = mydeal(num2cell(obj.regularity_r));
	[shp.regulari_x]  = mydeal(num2cell(obj.regularity_x));
	[shp.regulari_y]  = mydeal(num2cell(obj.regularity_x));
	[shp.runtime1]    = mydeal(num2cell(obj.runtime(:,1)));
	[shp.date]        = mydeal(obj.date);
	if (size(obj.runtime,2)>1)
		[shp.runtime2]   = mydeal(num2cell(obj.runtime(:,2)));
	end
	[shp.waveleng_r]  = mydeal(num2cell(obj.wavelength_r));
	[shp.waveleng_x]  = mydeal(num2cell(obj.wavelength_x));
	% note bug in shapewrite, logical variables are not exported
	[shp.qc_passed]  = mydeal(num2cell(double(obj.quality_check())));
	
	shp = Shp.set_geometry(shp,'Point');
	shp = Shp.nan2zero(shp);
	Shp.write(shp,shpname);
end % export_shp

