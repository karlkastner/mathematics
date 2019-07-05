% 2018-09-21 21:28:04.839819466 +0200
%% resample a track so that velocity vectors do not run into each other
function [x_,y_,u_,v_,S,S_] = resample_vector(x,y,u,v,dsmin,scale)
	method = 'linear';
	method = 'linear';
if (0)
	% iteratively
	x_ = x(1);
	y_ = y(1);
	dx = diff([x;x(1)]);
	dy = diff([y;y(1)]);
	ds = hypot(dx,dy);
	S  = cumsum(ds);
	S_ = S(1);

	idx=1;
	while (S_(end)<S(end))
		% get velocity
		u_(idx,1) = interp1(S,u,S_(idx),method);
		v_(idx,1) = interp1(S,v,S_(idx),method);
		% guess next point
		x_(idx+1,1) = interp1(S,x,S_(idx)+1,method);
		y_(idx+1,1) = interp1(S,y,S_(idx)+1,method);

		%while (true)
		dx = x_(idx+1,1) - x_(idx);
		dy = y_(idx+1,1) - y_(idx);

		% segment chord length
		ds = hypot(dx,dy);
		% projected arrow length
		la = scale*hypot(u_(idx)*dx/ds,v_(idx)*dy/ds);
		la = max(dsmin,la);
		% correct position
		dx = dx*la/ds;
		dy = dy*la/ds;
		x_(idx+1,1) = x_(idx)+dx;
		y_(idx+1,1) = y_(idx)+dy;
		S_(idx+1,1) = S_(idx)+la;
			% if segment length with arrow length tolerance
		%	if ( abs(ds-la) < 0.01;  )
		%		break;
		%	end
		%end
		idx=idx+1;
	end
	
	%u_(idx+1,1) = interp1(S,u,S_(idx));
	%v_(idx+1,1) = interp1(S,v,S_(idx));
	S_ = S_(1:end-1);
	x_ = x_(1:end-1);
	y_ = y_(1:end-1);
	S_ = S_*S(end)/S_(end);
else
	dx = diff([x;x(1)]);
	dy = diff([y;y(1)]);
	ds = hypot(dx,dy);
	S  = cumsum(ds);
	% direction
	dx = dx./ds;
	dy = dy./ds;
	% projected s of velocity
	% TODO, at segment centre
	umid =  u; %[mid(u); 0.5*(u(1)+u(end))];
	vmid =  v; %[mid(v); 0.5*(v(1)+v(end))];
	%vds = abs(dx.*umid + dy.*vmid);
%	vds = hypot(dx.*umid,dy.*vmid);
	vds = dx.*umid + dy.*vmid;

	S_ = 0;
	p  = 1 ;
	while (S_(end)<S(end))
		% TODO use interpolation object
		vds_ = interp1(S,vds,S_(end),method);
		S_(end+1,1) = S_(end) + max(dsmin,p*scale*abs(vds_));
		if (vds_ < 0)
			% implicit if pointing velocity is against
			for idx=1:10
			vds_ = interp1(S,vds,S_(end),method);
			S_(end,1) = S_(end-1) + max(dsmin,p*scale*abs(vds_));
			end
		end
		
		%S_(end+1,1) = S_(end) + max(dsmin,0.5*scale*(interp1(S,vds,S_end_) + interp1(S,vds,S_(end))));
	end
	% scale to identical length
	S_ = S_*S(end)/S_(end);
	x_ = interp1(S,x,S_,method);
	y_ = interp1(S,y,S_,method);
	u_ = interp1(S,u,S_,method);
	v_ = interp1(S,v,S_,method);
end
end % resample_vector

