% Di 12. Jan 11:02:39 CET 2016
%% Daniell window for smoothing the power spectrum
%% c.f. Daniell 1946
%% Bloomfield 2000
%% meko 2015
function f = daniell_window(m)
	if (m < 2)
		f = 1;
	else
		f = 1/(m-1)*ones(m,1);
		f(1) = 0.5*f(1);
		f(end) = 0.5*f(end);
	end
end

