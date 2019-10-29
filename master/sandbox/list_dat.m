% Tue May 29 16:17:09 MSK 2012
% Karl Kästner, Berlin

function list_dat(folder)
	if (nargin() < 1 || isempty(folder))
		folder = {'dat', 'dat-new'};
	end

	for fdx=1:length(folder)

		[flag name] = system(['ls -1 ' folder{fdx} '/fem-*.mat' ]);

		if (0 == flag)
			name = regexp(name, '\n', 'split');
			for idx=1:length(name)
				try
					s=load(name{idx});
					disp(sprintf('%s order %d L0 %s x0 %s k %d abstol %e', name{idx}, s.opt.order, ...
						mat2str(s.L0), mat2str(s.x0), s.k, s.opt.abstol ));
					clear s;
				catch
					disp([name{idx} ' error'])
				end
			end
		end % if 0 == flag

	end % for folder

	% FDM
	for fdx=1:length(folder)

		[flag name] = system(['ls -1 ' folder{fdx} '/fdm-*.mat' ]);

		if (0 == flag)
			name = regexp(name, '\n', 'split');
			for idx=1:length(name)
				try
					s=load(name{idx});
					disp(sprintf('%s dimension %d L0 %s x0 %s k %d abstol %e', ...
						name{idx}, s.d, mat2str(s.L0), mat2str(s.x0), s.k, s.abstol ));
					clear s;
				catch
					disp([name{idx} ' error'])
				end
			end
		end % if 0 == flag

	end % for folder
	
end % list_dat

