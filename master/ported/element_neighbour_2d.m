% Tue May  1 19:38:25 MSK 2012
% Karl Kästner, Berlin
%
% generate iverse mapping from triangle boundaries to triangles
function [N] = element_neighbour_2d(P, T, B, testflag)
	import java.util.Hashtable;
	
	if (nargin() < 4 || isempty(testflag))
		testflag = true;
	end

	% get length
	lt = size(T,1);
	lb = size(B,1);

	% allocate memory
	N = zeros(lt,3);

	% hash of triangle sides
	% key: lower_point_index, higher_point_index
	% value: triangle index, index of opposit point
	% if boundary, boundary index and negative number of boundary
	Sh = java.util.Hashtable(3*lt);

	% push the boundaries
	for idx=1:size(B,1)
		% push row and column in N, negative, domain
		Sh.put(hashkey(B(idx,1), B(idx,2)), -idx);
	end

	for idx=1:size(T,1)
		for jdx=1:3
			% point at opposit side
			key = hashkey(T(idx,mod(jdx,3)+1), T(idx,mod(jdx+1,3)+1));
			val = Sh.remove(key);
			if (isempty(val))
				Sh.put(key,[idx jdx]);
			else
				% Neighbourhood relation found
				%if (~isempty(val))
				N(idx,jdx) = val(1);
				if (val(1) > 0)
					N(val(1),val(2)) = idx;
				end
			end
		end % for jdx
	end % for idx
	if (testflag && ~Sh.isEmpty())
		K = Sh.keySet().toArray();
		for idx=1:length(K)
			disp( (Sh.get(K(idx)))');
		end
		error('fem_2d_element_boundary','inconsistent mesh');
	end
end % function fem_2d_element_boundary

