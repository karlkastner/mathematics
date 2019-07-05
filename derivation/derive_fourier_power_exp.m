% Tue  1 May 13:43:43 CEST 2018

function cp_C = derive()

syms c0 real
syms c1 c2
c=[c0 c1 c2]
o = [0 1 2];

% third power
for pdx=2:3
	cp = sym(zeros(2*pdx+1,1));

if (2 == pdx)
	for id1=1:length(c)
	for id2=1:length(c)
		[cp_12, o_12] = complex_exp_product([c(id1),c(id2)],[o(id1),o(id2)])
		for id5=1:2
			cp(o_12(id5)+1) = cp(o_12(id5)+1) + cp_12(id5);
		end
	end % id3
	end % id2
else
	for id1=1:length(c)
	for id2=1:length(c)
		[cp_12, o_12] = complex_exp_product([c(id1),c(id2)],[o(id1),o(id2)])
		for id3=1:length(c)
			for id4=1:2
				[cp_123, o_123] = complex_exp_product([cp_12(id4),c(id3)],[o_12(id4),o(id3)])
				for id5=1:2
					cp(o_123(id5)+1) = cp(o_123(id5)+1) + cp_123(id5);
				end
			end
		end % id3
	end % id2
end % id1
end
cp_C{pdx} = simplify(cp)

end % for pdx

cp_C{:}


%function step()
%	for id1=1:length(c)	
%	end
%end

end

