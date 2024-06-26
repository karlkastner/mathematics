% interpolation on fem
ipol:
	for 1:n in root
		pi = ipol_(T(i));
		if (null != pi) return
	end
end
ipol_:
	if (no_child)
		interpolate
		return 1
	end
	for each child
		s = get weight of linear combination
		% if max(s)
		if 1 == sum(s)
			return ipol(child)
		end
	end
	return 0
end

