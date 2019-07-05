% Thu 16 Feb 16:02:16 CET 2017
%% peaks of a periodogram
function [ma mi] = peaks_man(y,s)
	n = length(y);
	np = 0;
	id = [];
	lsig = [];
	rsig = [];
	for idx=2:n-1
		% check if this is a local peak
		if (  (1==idx || y(idx)>y(idx-1)) ...
                    & (n==idx || y(idx)>y(idx+1)))
		np = np+1;
		id(np) = idx;

		% find next max to the right
		[ldx] = find(y(1:idx-1) > y(idx),1,'last');
		if (isempty(ldx))
			ldx=0;
		end

		% find minimum in between (TODO minimum effect)
		[mv midx] = min(y(ldx+1:idx-1));
		midx = midx+ldx;
		
		% effect size
		ldy(np) = y(idx)-y(midx);

		% significance
		lsig(np) = abs(ldy(np)) - s(idx) - s(midx);
		%[y(idx), y(idx)-s(idx), y(midx), y(midx)+s(midx) lsig(np)]

		% find next max to the left
		rdx = idx + find(y(idx+1:end) > y(idx),1,'first');
		if (isempty(rdx))
			rdx=n+1;
		end
		
		% find minimum in between (TODO minimum effect)
		[mv midx] = min(y(idx+1:rdx-1));
		midx=midx+idx;
		
		% effect size
		rdy(np) = y(idx)-y(midx);

		% significance
		rsig(np) = abs(rdy(np)) - s(idx) - s(midx);
		end % if local max
	end % for idx
	ma.id  = id;
	ma.lsig = lsig;
	ma.rsig = rsig;

	% combine significance
	ma.sig = min(ma.lsig,ma.rsig);
	% sort
	[void sdx] = sort(ma.sig);
	%ma=rmfield(ma,'n');
	field_C = fieldnames(ma);
	for idx=1:length(field_C)
		ma.(field_C{idx}) = ma.(field_C{idx})(sdx);
	end % idx
end % peaks_man

