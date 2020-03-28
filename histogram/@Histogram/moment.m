% Di 8. Mär 09:38:23 CET 2016
% Karl Kastner, Berlin 
%
% non-central histogram moment
function [m obj] = moment(obj, order,varargin)
	m = obj.momentS(obj.h, obj.edge, order, varargin{:});
end

