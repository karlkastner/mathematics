% 2016-03-01 20:26:08.428033048 +0100
% Karl Kastner, Berlin
%
%
function [C2, level, midP] = quantile_envelope(X,Y,Z,x0,y0,p)
	n = length(x0);
	lowerZ = max(Z(:));
	lowerP = 0;
	upperZ = min(Z(:));
	upperP = 1;
	dz   = 0.5*abs(upperZ - lowerZ);
	midZ = 0.5*(upperZ + lowerZ);
	while (true)
		% compute mid-level contour
%		midZ = p_*upperZ+(1-p_)*lowerZ;
		C = contourc(X,Y,Z',[lowerZ,midZ,upperZ]);
		% extract centre-contour
		[Cc level] = extract_contour(C);
		% mid contour
		m = 0;
		C2 = {};
		% test which samples are contained
		for jdx=1:length(level)
			if (abs(level(jdx)-midZ) < 1e-7)
				% there can be several modes, count each of them
				C_   = [Cc{jdx}(1,:)',Cc{jdx}(2,:)'];
				flag = inpoly([x0,y0],C_);
				% count number of points contained in contour
				m = m+sum(flag);
				C2{end+1} = Cc{jdx};
			end
		end
		midP = m/n;
		dz = dz/2;
%		[midZ dz midP]
%pause
		%[lowerZ midZ upperZ]
		if (midP < p)
			midZ = midZ - dz;
		else
			midZ = midZ + dz;
		end

%		[lowerP midP upperP]

		% stop if sufficiently close
		% TODO, this may fail if points have same distance
		if ( abs(midP-p) < 1.5/n)
			break;
		end
		% discard furthest end
		if ( upperP-p > p-lowerP )
			% do not move the upper P, but just shift midz
			upperP = midP;
			upperZ = midZ;
		else
			lowerP = midP;
			lowerZ = midZ;
		end
	end
end

