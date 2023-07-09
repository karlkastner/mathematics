% Thu  6 Jul 09:32:32 CEST 2023
function init(obj)
	obj.stat.area_msk = NaN;
	obj.stat.contrast = NaN;
	obj.stat.coverage = NaN;
	obj.stat.isisotropic = NaN;
	obj.stat.stati.intS_hp_sig = NaN;
	field_C = {'x','y','r'};
	for field = field_C
		obj.stat.L_eff.(field{1}) = NaN;
	end
	obj.stat.p_periodic = NaN;
	obj.stat.p_S_hp = NaN;
	field_C = {'hat','hp','bar'};
	for field = field_C
		obj.stat.fc.radial.(field{1}) = NaN;
		obj.stat.fc.x.(field{1}) = NaN;
		obj.stat.Sc.radial.(field{1}) = NaN;
		obj.stat.Sc.x.(field{1}) = NaN;
		obj.stat.Sc.y.(field{1}) = NaN;
	end
end

