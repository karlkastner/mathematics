
% Mi 3. Feb 12:46:57 CET 2016
	function obj = init_measure_map(obj)
		nt_sample         = max(0,floor((obj.T(end)-obj.T(1))/obj.out.map.dt)-1);
		obj.out.map.dt  = (obj.T(end)-obj.T(1))/nt_sample;
		% allocate memory
		obj.out.map.A   = zeros(obj.nx,nt_sample);
		obj.out.map.Q   = zeros(obj.nx,nt_sample);
		obj.out.map.T   = zeros(nt_sample,1);
		% start counter
		obj.out.map.t   = obj.T(1);
		obj.out.map.idx = 0;
	end % init_measure_map
	function obj = measure_map(obj,t,A,Q)
		if (t >= obj.out.map.t+obj.out.map.dt)
			% index counter
			obj.out.map.idx                = obj.out.map.idx + 1;
			obj.out.map.T(obj.out.map.idx) = t;
			obj.out.map.A(:,obj.out.map.idx) = A;
			obj.out.map.Q(:,obj.out.map.idx) = Q;
			% time of next sample
			obj.out.map.t        = obj.out.map.t+obj.out.map.dt;
		end % if
	end % measure_map

	function [obj]   = init_measure_points(obj)
		nt_sample         = max(0,floor((obj.T(end)-obj.T(1))/obj.out.point.dt)-1);
		obj.out.point.dt  = (obj.T(end)-obj.T(1))/nt_sample;
		% determine output indices
		n_sample          = length(obj.out.point.x);
		obj.out.point.id  = interp1(obj.x,1:length(obj.x),obj.out.point.x,'nearest');
		obj.out.point.id  = min(length(obj.x),max(1,obj.out.point.id));
		% allocate memory
		obj.out.point.A   = zeros(n_sample,nt_sample);
		obj.out.point.Q   = zeros(n_sample,nt_sample);
		obj.out.point.T   = zeros(nt_sample,1);
		% start counter
		obj.out.point.t   = obj.T(1);
		obj.out.point.idx = 0;
	end % init_measure_points
	function [obj] = measure_points(obj,t,A,Q)
		if (t >= obj.out.point.t+obj.out.point.dt)
			% index counter
			obj.out.point.idx      = obj.out.point.idx + 1;
			obj.out.point.T(obj.out.point.idx) = t;
			obj.out.point.A(:,obj.out.point.idx) = A(obj.out.point.id);
			obj.out.point.Q(:,obj.out.point.idx) = Q(obj.out.point.id);
			% time of next sample
			obj.out.point.t        = obj.out.point.t+obj.out.point.dt;
		end 
	end % measure_points
