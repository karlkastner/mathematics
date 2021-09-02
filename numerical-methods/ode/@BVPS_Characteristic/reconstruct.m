% Wed 26 Aug 16:27:45 +08 2020
% TODO reconstruct with exp instead of using inner2outer
function reconstruct(obj,ypm)
	nci = obj.nci;
	npii = obj.npii;
	ni = obj.ni;
	
	% TODO
	% yc = obj.Ic*ypm;

    for cdx=1:obj.nc
        yc = zeros(nci(end,cdx)-1,1);
	if (obj.opt.reconstruct_y)
		y = zeros(ni(end,cdx)-1,1);
	end
	for edx=1:obj.neq
		c = obj.out(cdx).cc;
		switch (obj.oo(edx))
		case {1}
		    if (~obj.opt.dischargeisvariable)
			% TODO, check for each row
		     if (0 ~= c(1,2,edx))
			% this does not happen for the
			yc_ = (   ypm(npii(edx,cdx)  :2:npii(edx+1,cdx)-2) ...
			       + ypm(npii(edx,cdx)+1:2:npii(edx+1,cdx)-1) ...
			     );
		     else
			% degenerated, linear function
			yc_ = ypm(npii(edx,cdx)+1:2:npii(edx+1,cdx)-1); % or +0, -2?
		     end
			if (obj.opt.reconstruct_y)
				y_ = inner2outer(yc_);
			end
		    else
		     if (0 ~= c(1,2,edx))
			yc_ = (  ypm(npii(edx,cdx)  :2:npii(edx+1,cdx)-3) ...
			       + ypm(npii(edx,cdx)+1:2:npii(edx+1,cdx)-2) ...
			      );
		     else
			% degenerated, linear function
			yc_ = ypm(npii(edx,cdx)+1:2:npii(edx+1,cdx)-2); % or +0, -2?
		     end
			if (obj.opt.reconstruct_y)
				y_  = [inner2outer(yc_); ypm(npii(edx+1,cdx)-1)];
			end
			yc_ = [yc_; ypm(npii(edx+1,cdx)-1)];
		    end
		case {2}
			% TODO use bvp1c/2c to expand with exact exponential
			yc_ = (  ypm(npii(edx,cdx)  :3:npii(edx+1,cdx)-3) ...
			       + ypm(npii(edx,cdx)+1:3:npii(edx+1,cdx)-2) ...
			       + ypm(npii(edx,cdx)+2:3:npii(edx+1,cdx)-1) ...
		              );
			if (obj.opt.reconstruct_y)
				y_ = inner2outer(yc_);
			end
	        otherwise
			error('here');
	     	end % switch mm
		yc(nci(edx,cdx):nci(edx+1,cdx)-1) = yc_;
		if (obj.opt.reconstruct_y)
			y(ni(edx,cdx):ni(edx+1,cdx)-1) = y_;
		end
	end % for edx
    	obj.out(cdx).yc   = yc;
	obj.out(cdx).y    = y;
	obj.out(cdx).ypm  = ypm(npii(1,cdx):npii(end,cdx)-1);
    end % for cdx
end % reconstruct

