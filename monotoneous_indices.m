% Tue 31 Jan 15:27:52 CET 2023
function fdx = monotoneous_indices(x,mode)
	x= cvec(x);
	switch (mode)
	case {'descending','d',-1}
		xmin = cummin(x);
		fdx  = [true;x(2:end)<xmin(1:end-1)];
	case {'ascending','a',1}
		xmax = cummax(x);
		fdx  = [true;x(2:end)>xmax(1:end-1)];
	otherwise
		error('unknown mode')
		disp('mode');
	end
end

