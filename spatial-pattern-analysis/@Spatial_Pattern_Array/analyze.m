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
		disp(sprintf('%s %d/%d Minutes spend %0.0f remaining %0.0f',obj.type,idx,n,t/60,t/(idx-1)*(n-idx+1)/60));
		dx = obj.dx_sample(idx);
		try
			if (~isfinite(dx)||(dx<=0))
				error('sampling resolution not specified');	
			end
			dx = min(obj.opt.dx_max,max(obj.opt.dx_min,dx));
			% generate filename
			[folder,filename,folder_a] = obj.generate_filename(idx,dx);
			mkdir([folder,filesep,'analysis']);
			spname = [folder_a,filesep,filename(1:end-4),'.mat'];
			obj.filename_C{idx} = filename;
			% check if image has already been analyzed
			if (obj.opt.reload && exist(spname,'file'))
				% reload
				clear sp
				load(spname,'sp');
				if (~isfield(sp.stat,'stati'))
					sp.imread([folder,filesep,filename]);
					sp.opt.test_for_periodicity = obj.opt.test_for_periodicity;
					sp.prepare_analysis();
					sp.analyze_grid();

					% clear 2d varibales to save space
					sp.clear_2d_properties();	
		
					%sp = sp_a(idx);
					mkdir(dirname(spname));
					runtime = toc(timer)-t;
					obj.runtime(idx,2) = runtime;
					% save analysis of individual pattern
					save(spname,'sp','runtime');
				end
			else
				% analyze pattern
				disp([num2str(idx), ' ', num2str(dx), ' ', filename])

				sp = Spatial_Pattern();
	
				sp.opt.test_for_periodicity = obj.opt.test_for_periodicity;
		
				sp.imread([folder,filesep,filename]);
		
				sp.analyze_grid();
		
				% clear 2d varibales to save space
				sp.clear_2d_properties();	
		
				%sp = sp_a(idx);
				mkdir(dirname(spname));
				runtime = toc(timer)-t;
				obj.runtime(idx,2) = runtime;
				% save analysis of individual pattern
				save(spname,'sp','runtime');

			end % else of (reload && exist spname)
			obj.sp_a(idx) = sp;
		catch e
			%disp(base);
			disp(e)
			s = e.stack;
			for idx=1:length(s)
				disp(s(idx));
			end	
			obj.error_C{idx,2} = e;
		end % of try-catch 
	end % for idx (each pattern)
end % analyze

