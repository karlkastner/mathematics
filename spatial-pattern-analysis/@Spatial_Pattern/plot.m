% Wed  8 Feb 09:58:38 CET 2023

function plot(obj,field_str)
	field_C  = strsplit(field_str,'.');
	var      = getfield_deep(obj,field_str);
	lambda_c = obj.lambda_c();
	if (~isfinite(lambda_c))
		lambda_c = 1./obj.stat.f_50;
	else
		lambda_c = lambda_c;
	end

	if (~isvector(var))
	switch (field_C{end})
	case {'b'}
		%s = size(img);
		%x = (0:obj.stat.siz(1));
		%y = (0:obj.stat.siz(2));
		if (isempty(obj.msk.b))
			imagesc((obj.x-obj.x(1))/lambda_c,(obj.y-obj.y(1))/lambda_c,obj.b_square);
		else
			%surface((obj.x-obj.x(1))/lambda_c,(obj.y-obj.y(1))/lambda_c,obj.b_square,'Edgecolor','none');
			imagesc((obj.x-obj.x(1))/lambda_c,(obj.y-obj.y(1))/lambda_c,obj.b_square,'AlphaData',obj.msk.b_square)
		end % else of isempty alpha
		%title(title_C{:});
		axis xy
		axis equal
		axis tight
		xlabel('$x/\lambda_c$','interpreter','latex');
		ylabel('$y/\lambda_c$','interpreter','latex');
	case {'hat','clip','bar'}
	switch (field_C{1})
	case {'S'}
		imagesc(fftshift(obj.f.x*lambda_c),fftshift(obj.f.y*lambda_c),fftshift(var'/lambda_c.^2));
		%axis(4*[-1,1,-1,1]);
		axis xy;
		s = 4;
		axis equal;
		axis tight
		xlim(s*[-1,1]);
		ylim(s*[-1,1]);
		colormap(flipud(gray));
		xlabel('Wavenumber $k_x/k_c$','interpreter','latex');
		ylabel('Wavenumber $k_y/k_c$','interpreter','latex');
	case {'R'}
	%if (pline)
	%	surface(fc(kdx)*obj.y,fc(kdx)*obj.x,obj.R.clip,'edgecolor','none')
	%else
		imagesc((obj.x)/lambda_c,(obj.y)/lambda_c,fftshift(var)');
	%end
	%hold on
	%if (lineflag)
		%plot(fc(kdx)*y_(round(jc([1,end]))),fc(kdx)*x_(round(ic([1,end]))),'r','linewidth',1.5);
	%end
		axis xy 
		axis equal
		axis tight
		axis([obj.opt.xlim,obj.opt.xlim]); % rp
	otherwise
		error(['Unknown field ', field_str]);
	end

	otherwise
		error(['Unknown field ', field_str]);
	end
	else
	switch (field_C{end-1})
	case {'radial'}
		plot(obj.f.r*lambda_c,var/lambda_c);
		xlabel('Wavenumber $k_r/k_c$','interpreter','latex');
		ylabel('Density $S_r/\lambda_c$','interpreter','latex');
		% TODO no magic numbers
		xlim(obj.opt.xlim);
	case {'x'}
		fdx = obj.f.x>=0;
		plot(obj.f.x(fdx)*lambda_c,var(fdx)/lambda_c);
		xlabel('Wavenumber $k_x/k_c$','interpreter','latex');
		ylabel('Density $S_x/\lambda_c$','interpreter','latex');
		xlim(obj.opt.xlim);
	case {'y'}
		fdx = obj.f.y>=0;
		plot(obj.f.y(fdx)*lambda_c,var(fdx)/lambda_c);
		xlabel('Wavenumber $k_y/k_c$','interpreter','latex');
		ylabel('Denisty $S_y / \lambda_c$','interpreter','latex');
		xlim(obj.opt.xlim);
	case {'angular'}
		fdx = obj.f.angle >= -pi/2 & obj.f.angle <= pi/2;
		plot(obj.f.angle(fdx),var(fdx));
		xlabel('Angle $\theta / \pi$','interpreter','latex');
		ylabel('Density $S_\theta$','interpreter','latex');
		set(gca,'xtick',[0,1/4,1/2]*pi,'xticklabel',{'0','\pi/4','\pi/2'});
		%xlim([-pi,pi]/2);
		xlim([-1,1]*pi/2);
	otherwise
		error(['Unknown field ', field_str]);
	end % switch field_C{end-1}
	end % of else ~isvector
end %

