% 2014-06-29 15:39:10.090386438 +0200
% Karl Kastner, Berlin
%
%% get first intersection between lines in A and B
function [flag, idx_, jdx_] = first_intersect(A,B)
	flag = 0;
	idx_ = [];
	jdx_ = [];
	for idx=1:size(A,1)-1
	 for jdx=1:size(B,1)-1
		if (intersect(A(idx,:),A(idx+1,:),B(jdx,:),B(jdx+1,:)))
			idx_ = idx;
			jdx_ = jdx;
			return;
		end
         end
	end
	% TODO, wrap first and last elements
