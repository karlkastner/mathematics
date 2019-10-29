
% until HH
'HH'
m = 0;
n = 1000000;
for idx=1:n
	T_old = randn();
	while (1)
		m = m+1;		
		T = randn();
		if (T_old <= 0 && T <= 0)
			break;
		end
		T_old = T;
	end
end

m/n

'TT'
m = 0;
n = 1000000;
for idx=1:n
	T_old = randn();
	while (1)
		m = m+1;		
		T = randn();
		if (T_old >= 0 && T >= 0)
			break;
		end
		T_old = T;
	end
end

m/n


% until TH
'TH'

m = 0;
for idx=1:n
	T_old = randn();
	while (1)
		m = m+1;		
		T = randn();
		if (T_old <= 0 && T >= 0)
			break;
		end
		T_old = T;
	end
end

m/n

% until HT
'HT'

m = 0;
for idx=1:n
	T_old = randn();
	while (1)
		m = m+1;		
		T = randn();
		if (T_old <= 0 && T >= 0)
			break;
		end
		T_old = T;
	end
end

m/n

% how often in sequence
'seq'

'HT'

m = 0;
o=0;
T_old = randn();
for idx=1:n
	while (1)
		T = randn();
		if (T_old <= 0 && T >= 0)
			m = m+1;		
			break;
		end
		if (T_old <= 0 && T <= 0)
			o = o+1;		
			break;
		end
		T_old = T;
	end
end

m/n
o/n

