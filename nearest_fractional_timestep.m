% Tue 24 Nov 15:28:25 +08 2020
function dt = nearest_fractional_timestep(dt,T)
	for idx=[1,2:2:100]
		dt_ = T/idx;
		if (dt_ < dt)
			dt = dt_;
			break;
		end
	end
end

