% 2016-09-22 19:55:22.401333352 +0200 nl.m
% 2016-08-12 10:32:20.605122006 +0200
% Karl Kastner, Berlin
%
%% number rows (lines) of a matrix
%%
%% analogue to unix nl command
function x = nl(x)
	if (~iscell(x))
		x = [(1:size(x,1))',x];
	else
		x(:,2:end+1) = x;
		x(:,1) = num2cell((1:size(x,1))');
	end
end

