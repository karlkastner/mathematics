% 2023-07-26 13:43:55.873240646 +0200
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
%% determine the sampling interval for fetching images from the Google satellite
%% server and later processing
%%
function obj = fetch(obj,inshpname,outshpname)
    shp   = Shp.read(inshpname);
    % optionally skip patterns for quick testing
    shp = shp(1:obj.opt.skip:end);
    % important, otherwise masks are not burned
    shp   = Shp.reassign_id(shp);
    n     = length(shp);
    obj.centroid = Shp.centroid(shp);

    % initialize with dx0
    % obj.dx_sample       = repmat(obj.opt.dx0,obj.n,1);
    [shp.dx_sample]     = deal(obj.opt.dx0);
    [shp.fr_max_target] = deal(-1);
    [shp.filename]      = deal('');
    % write as data is passed this way to the python script
    Shp.write(shp,outshpname);

    level = 1;
    obj.level_a = ones(obj.n,1);
    % dx = obj.opt.dx0;

    % iteratively decrease resolution until lowest 
    timer = tic();
    %while (dx >= obj.opt.dx_min)
    while (true) % obj.resolution(level) >= obj.opt.dx_min)
        % obj.resolution(level) = dx;
    	dx = obj.resolution(level);
    
        % fetch patterns that are not yet well resolved at current resolution
        cmd_str = sprintf(obj.opt.cmd_str,obj.type,dx);
        system(cmd_str);
    
        % for each pattern
        for idx=1:obj.n
            t0 = toc(timer);
            % only process, if not yet fully resolved at the current resolution level
            if (obj.level_a(idx) >= level)
            %if (obj.dx_sample(idx) <= dx)
    
                [folder,filename,folder_a] = obj.generate_filename(idx,level);
                obj.filename_C{idx}    = filename;
                %mkdir([dirname(folder),filesep,'analysis',filesep,basename(folder)]);
    
                spname = [folder_a,filesep,filename(1:end-4),'.mat'];
                disp([obj.type, ' ', num2str(dx), ' ', num2str(idx), ' quantile estimation ', filename]);
                
                try
                    if (obj.opt.reload && exist(spname))
                        load(spname,'sp');
                    else
                        % check that file is below size limit
                        info=imfinfo([folder,filesep,filename]);
                        A = max(info.Width,info.Height)^2;
                        if (A > obj.opt.area_max)
                            error('Image exceeds size limit');
                        end
        
                        sp = Spatial_Pattern();
                        sp.imread([folder,filesep,filename]);
                        sp.prepare_analysis();
                        % clear 2d varibales to save space
                        sp.clear_2d_properties();
                        runtime = toc(timer)-t0;
                        % note, that only stores runtime of pattern preparation in highest resolution
                        obj.runtime(idx,1) = runtime;
                        mkdir(folder_a);
                        save(spname,'sp');
                    end % ~isfield
                    % TODO : use simpler requirement, like q90 < frmax/2, tail is anyway not normal
                    obj.fr_max_target(idx,1) = (sp.stat.q.fr.p50 + obj.opt.tail_trim_scale*(sp.stat.q.fr.p84-sp.stat.q.fr.p50));
                    %disp([idx,dx,s.(dx_field).q.fr.p50, s.(dx_field).q.fr.p84, s.(dx_field).q.fr.p95, s.(dx_field).q.fr.max, s.fr_max_target < s.(dx_field).q.fr.max]);
                    if (obj.fr_max_target(idx) > sp.stat.q.fr.max);
                        % undersampled
                        % obj.dx_sample(idx,1) = dx/2;
			obj.level_a(idx,1) = level+1;
                    elseif (obj.fr_max_target(idx) < 0.5*sp.stat.q.fr.max);
                        % oversampled
                        % obj.dx_sample(idx,1) = 2*dx;
			obj.level_a(idx,1) = level-1;
                    else
                        % good sampling range
                        % obj.dx_sample(idx,1) = dx;
			obj.level_a(idx,1) = level;
                    end
                    shp(idx).dx_sample = obj.resolution(obj.level_a(idx));
		    % obj.dx_sample(idx);
                catch e
                    disp(e);
                    for sdx=1:length(e.stack)
                        disp(e.stack(sdx))
                    end
                    disp(obj.centroid(idx,:));
                    %obj.dx_sample(idx) = NaN;
		    obj.level_a(idx) = NaN;
                    shp(idx).dx_sample = -1;
                    obj.error_C{idx,1} = e;
                end % catch of try
            % this must be at the end of the catch block, in case something goes wrong
            end % if dx_sample <= dx
        end % for idx (each pattern)
        Shp.write(shp,outshpname);
        % dx    = dx / 2;
        level = level+1;
	if (obj.resolution(level) < obj.opt.dx_min)
		break;
	end
    end % while dx >= dx_min
    % final fetching for coarse resolution patterns
    if (obj.opt.dx_min < obj.opt.dx_max)
        cmd_str = sprintf(obj.opt.cmd_str,obj.type,obj.opt.dx_max);
        system(cmd_str);
    end
end % fetch

