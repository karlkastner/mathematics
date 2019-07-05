% 2015-06-26 18:59:35.306986477 +0200
function varargout = transposeall(func,varargin)
	C=cell(nargout(func),1);
	[C{:}] = func(varargin{:});
%	varargin{:}
	varargout = cellfun(@transpose,C,'uniformoutput',false);
%	varargout = cell(size(C))
%	C
%	[varargout{:}] = deal(C{:});
%	 [a b] = deal(C{:})
%{mean_man([1:10; (1:10).^2])},
%	for idx=1:length(C)%varargin)
%		varargout{idx} = varargin{idx}';
%	end
end

