% Wed 19 May 13:30:08 CEST 2021
function [To,y] = euler(fun,To,y0,dt)
	nt = round((To(end)-To(1))/dt);
	t  = linspace(To(1),To(end),nt+1);

	if (length(To)<=2)
		To = t;
	end
	nto = length(To);
	y   = zeros(length(y0),nto);
	y(:,1) = y0;
	odx = 1;
	for idx=2:nt+1
		if (t(idx)>=To(odx))
			y(:,odx) = y0;
			odx=odx+1;
		end
		dt    = (t(idx)-t(idx-1));
		dy_dt = fun(t,y0);
		y0 = y0 + dt*dy_dt;
	end
	y(:,end) = y0;
	y = y.';
end

