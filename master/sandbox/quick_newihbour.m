% Wed Jun 20 15:50:44 MSK 2012

% neighbourhood reference
lt = size(T,1);
I = (1:lt)';
edge_buf = [ [min([T(:,1) T(:,2)],[],2) max([T(:,1) T(:,2)],[],2)] I;
	  [min([T(:,1) T(:,3)],[],2) max([T(:,1) T(:,3)],[],2)] I;
	  [min([T(:,2) T(:,3)],[],2) max([T(:,2) T(:,3)],[],2)] I];
edge_mat = sparse(edge_buf);
N = [   [edge_mat(size(edge_mat, sub2ind(min([T(:,1) T(:,2)],[],2), max([T(:,1) T(:,2)],[],2))) - I]
	[edge_mat(size(edge_mat, sub2ind(min([T(:,1) T(:,3)],[],2), max([T(:,1) T(:,3)],[],2))) - I]
	[edge_mat(size(edge_mat, sub2ind(min([T(:,2) T(:,3)],[],2), max([T(:,2) T(:,3)],[],2))) - I] ];

