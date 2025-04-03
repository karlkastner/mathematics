% Wed  8 Feb 09:58:38 CET 2023
% Karl Kästner, Berlin
%
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
%
%% plot the pattern, its density or autocorrelation
function [c,cbh] = plot(obj,field_str,varargin)
	field_C  = strsplit(field_str,'.');
	var      = getfield_deep(obj,field_str);
	lambda_c = obj.lambda_c();
	if (~isfinite(lambda_c))
		lambda_c = 1./obj.stat.q.fr.p50;
	else
		lambda_c = lambda_c;
	end
	if (~isvector(var))
	switch (field_C{1})
	case {'b'}
		if (isempty(obj.msk.b))
			imagesc((obj.x-obj.x(1))/lambda_c,(obj.y-obj.y(1))/lambda_c,obj.b_.square);
		else
			%surface((obj.x-obj.x(1))/lambda_c,(obj.y-obj.y(1))/lambda_c,obj.b_square,'Edgecolor','none');
			imagesc((obj.x-obj.x(1))/lambda_c,(obj.y-obj.y(1))/lambda_c,obj.b_.square,'AlphaData',obj.msk.b_square)
		end % else of isempty alpha
		axis xy
		axis equal
		axis tight
		xlabel('Position $x/\lambda_c$','interpreter','latex');
		ylabel('Position $y/\lambda_c$','interpreter','latex');
	case {'msk'}
	switch (field_C{2})
	case {'f'}
		imagesc(fftshift(obj.f.x*lambda_c),fftshift(obj.f.y*lambda_c),fftshift(obj.msk.f));
		axis xy;
		axis equal;
		axis tight
		s = 10;
		xlim(s*[-1,1]);
		ylim(s*[-1,1]);
		colormap(flipud(gray));
		xlabel('Wavenumber $k_x/k_c$','interpreter','latex');
		ylabel('Wavenumber $k_y/k_c$','interpreter','latex');
	end
	case {'S'}
		imagesc(fftshift(obj.f.x*lambda_c),fftshift(obj.f.y*lambda_c),fftshift(var/lambda_c.^2));
		axis xy;
		axis equal;
		axis tight
		s = 4;
		xlim(s*[-1,1]);
		ylim(s*[-1,1]);
		colormap(flipud(gray));
		xlabel('Wavenumber $k_x/k_c$','interpreter','latex');
		ylabel('Wavenumber $k_y/k_c$','interpreter','latex');
		c = colorbar();
		if (strcmp(field_str,'S.hat') || strcmp(field_str,'S.rot.hat') || strcmp(field_str,'S.rot.hp'))
			title(c,'$\hat S/\lambda_c^2$','interpreter','latex');
		else
			title(c,'$S/\lambda_c^2$','interpreter','latex');
		end
	case {'R'}
	%if (pline)
	%	surface(fc(kdx)*obj.y,fc(kdx)*obj.x,obj.R.clip,'edgecolor','none')
	%else
		imagesc((obj.x)/lambda_c,(obj.y)/lambda_c,fftshift(var));
	%end
	%hold on
	%if (lineflag)
		%plot(fc(kdx)*y_(round(jc([1,end]))),fc(kdx)*x_(round(ic([1,end]))),'r','linewidth',1.5,varargin{:});
	%end
		xlabel('Lag distance $x/\lambda_c$','interpreter','latex');
		ylabel('Lag distance $y/\lambda_c$','interpreter','latex');
		axis xy 
		axis equal
		axis tight
		axis([obj.opt.xlim,obj.opt.xlim]); % rp
		c = colorbar();
		if (strcmp(field_str,'R.hat'))
			title(c,'$\hat R$','interpreter','latex');
		else
			title(c,'$R$','interpreter','latex');
		end
%	otherwise
%		error(['Unknown field ', field_str]);
%	end

	otherwise
		error('Spatial_Pattern:UnknownField', ['Unknown field ', field_str]);
	end
	else
        switch (field_C{1})
	case {'S'}
	switch (field_C{end-1})
	case {'radial'}
		plot(obj.f.r*lambda_c,var/lambda_c,varargin{:});
		xlabel('Wavenumber $k_r/k_c$','interpreter','latex');
		ylabel('Density $S_r/\lambda_c$','interpreter','latex');
		% TODO no magic numbers
		xlim(obj.opt.xlim);
	case {'x'}
		plot(fftshift(obj.f.x)*lambda_c,fftshift(var)/lambda_c,varargin{:});
		xlabel('Wavenumber $k_x/k_c$','interpreter','latex');
		ylabel('Density $S_x/\lambda_c$','interpreter','latex');
		xlim(obj.opt.xlim);
	case {'xp'}
		plot(fftshift(obj.f.x)*lambda_c,fftshift(var)/lambda_c,varargin{:});
		xlabel('Wavenumber $k_x/k_c$','interpreter','latex');
		ylabel('Density $S_x^+/\lambda_c$','interpreter','latex');
		xlim(obj.opt.xlim);
	case {'y'}
		plot(fftshift(obj.f.y)*lambda_c,fftshift(var)/lambda_c,varargin{:});
		xlabel('Wavenumber $k_y/k_c$','interpreter','latex');
		ylabel('Denisty $S_y / \lambda_c$','interpreter','latex');
		xlim(obj.opt.xlim);
	case {'angular'}
		plot(obj.f.angle,var,varargin{:});
		xlabel('Angle $\theta$','interpreter','latex');
		ylabel('Density $S_\theta$','interpreter','latex');
		set(gca,'xtick',[-1/2,-1/4,0,1/4,1/2]*pi,'xticklabel',{'-\pi/2','-\pi/4','0','\pi/4','\pi/2'});
		%xlim([-pi,pi]/2);
		xlim([-1,1]*pi/2);
		ylim([0,ceil(10.1*max(var))/10]);
	case {'angular_p'}
		fdx = obj.f.angle >= -pi/2 & obj.f.angle <= pi/2;
		plot(obj.f.angle(fdx),var(fdx),varargin{:});
		xlabel('Angle $\theta$','interpreter','latex');
		ylabel('Density $S_\theta^+$','interpreter','latex');
		set(gca,'xtick',[-1/2,-1/4,0,1/4,1/2]*pi,'xticklabel',{'-\pi/2','-\pi/4','0','\pi/4','\pi/2'});
		%xlim([-pi,pi]/2);
		xlim([-1,1]*pi/2);
		ylim([0,ceil(10.1*max(var(fdx)))/10]);
	otherwise
		error('Spatial_Pattern:UnknownField',['Unknown field ', field_str]);
	end % switch field_C{end-1}
	case {'R'}
	switch (field_C{end-1})
	case {'radial'}
		plot(obj.r/lambda_c,var,varargin{:});
		xlabel('Lag distance $x/\lambda_c$','interpreter','latex');
		ylabel('Autocorrelation $R_r$','interpreter','latex');
		% TODO no magic numbers
		xlim(obj.opt.xlim);
	otherwise
		error('Spatial_Pattern:UnknownField',['Unknown field ', field_str]);
	end % switch field_C{end-1}
	end
	end % of else ~isvector
end %

