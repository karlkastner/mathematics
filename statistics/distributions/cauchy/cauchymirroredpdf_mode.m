% 2024-06-30 11:55:28.738845211 +0200
function [fc,Sc] = cauchymirroredpdf_mode(f0,s)
	   r = 2*f0*sqrt(f0^2 + s^2) - f0^2 - s^2
	   if (r >= 0)
	   	fc = sqrt(r);
	   else
		r = - 2*f0*sqrt(f0^2 + s^2) - f0^2 - s^2
		if (r >= 0)
		   	fc = sqrt(r);
		else
			fc = 0;
		end
	   end
	Sc = cauchymirroredpdf(fc,f0,s);
end

