% Di 2. Feb 11:07:09 CET 2016
% Karl Kastner, Berlin
%
%% effective degrees of freedom for weighted samples
% see: Kish, survey sampling
% TODO, this is identical to effectife_sample_size
function dof = wdof(w)
	dof = sum(w).^2/sum(w.^2);
end

