% Wed 15 Feb 18:00:11 CET 2023
function S2d = reconstruct_isotropic_density(fr,angle,Sr,Sa)
	if (isa(Sr,'function_handle'))
		Sr = Sr(fr);
	end
	if (isa(Sa,'function_handle'))
		Sa = Sa(angle);
	end
	S2d = 0.25*Sr.*Sa;
	S2d(~isfinite(S2d)) = 0;
end

