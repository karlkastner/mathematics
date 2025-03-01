% 2022-09-26 14:32:10.449630621 +0200
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
%% analyze spatial patterns
function analyze(obj)
	n    = size(obj.centroid,1);
	timer = tic();

	% process each pattern
	for idx=1:n
 		t = toc(timer);
		str = sprintf('%s Pattern %d/%d Time %0.0f min spend %0.0f min remaining',obj.type,idx,n,t/60,t/(idx-1)*(n-idx+1)/60);
		% dx    = obj.dx_sample(idx);
		level = obj.level_a(idx);
		try
			clear sp
			if (~isfinite(level)) % isfinite(dx) || (dx<=0))
				error('Spatial_Pattern_Array:SamplingResolutionNotSpecified','Sampling resolution not specified');	
			end
			level = min(max(level,obj.level_min),obj.level_max);
			% dx = min(obj.opt.dx_max,max(obj.opt.dx_min,dx));
			% the next line should not be necessary
			% TODO apply to level
			% generate filename
			[folder,filename,folder_a] = obj.generate_filename(idx,level);
			mkdir(folder_a); %[folder,filesep,'analysis']);
			spname = [folder_a,filesep,filename(1:end-4),'.mat'];
			obj.filename_C{idx} = filename;
			% check if image has already been analyzed
			if (obj.opt.reload && exist(spname,'file'))
				% reload
				clear sp
				load(spname,'sp');
				% analyze pattern if output file exists but
				% pattern not yet analyzed, this is the case 
				% after determining the sampling resolution
				% at the first run
				if (~sp.isanalyzed());
					dx    = obj.resolution(level);
					disp([str,' analyzing ',num2str(dx)]);
					analyze_();
				else
					disp([str,' reloading']);
				end
			else
				% analyze pattern
				dx    = obj.resolution(level);
				sp = Spatial_Pattern();
				sp.init();
				if (obj.opt.analyze)
					disp([str,' analyzing ',num2str(dx)]);
					analyze_();
				else
					disp([str, ' skipping']);
				end
			end % else of (reload && exist spname)
		catch e
			sp = Spatial_Pattern();
			sp.init();
			disp(e);
			s = e.stack;
			for jdx=1:length(s)
				disp(s(jdx));
			end	
			obj.error_C{idx,2} = e;
		end % of try-catch
		obj.sp_a(idx,1) = sp;
	end % for idx (each pattern)

function analyze_()
		info=imfinfo([folder,filesep,filename]);
		A = max(info.Width,info.Height)^2;
		if (A > obj.opt.area_max)
			error('Spatial_Pattern_Array:ImageExceedsSizeLimit','Image exceeds size limit');
		end

		sp.imread([folder,filesep,filename]);
		sp.opt.test_for_periodicity = obj.opt.test_for_periodicity;

		sp.analyze_grid();
		try
			sp.fit_parametric_densities();
		catch e
			disp(e);
			disp(e.identifier);
			for edx=1:length(e.stack);
				disp(e.stack(edx));
			end
			obj.error_C{idx,3} = e;
		end

		% clear 2d varibales to save space
		sp.clear_2d_properties();
		sp.cast(@single);	

		mkdir(dirname(spname));

		% save analysis of individual pattern
		save(spname,'sp');
end

end % analyze

