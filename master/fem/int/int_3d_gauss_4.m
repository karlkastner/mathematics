% Wed Jul 11 18:11:58 MSK 2012
% Wed Jul 11 17:49:26 MSK 2012
% Karl Kästner, Berlin

% p = 4 (fehler hier, sicher p nur 2)
function [w b flag] = int_3d_gauss_4()
	w = [0.250000000000000; 0.250000000000000; 0.250000000000000; 0.250000000000000];
	b = [   0.138196601125011   0.138196601125011   0.138196601125011   0.585410196624969;
		0.138196601125011   0.138196601125011   0.585410196624969   0.138196601125011;
		0.138196601125011   0.585410196624969   0.138196601125011   0.138196601125011;
		0.585410196624969   0.138196601125011   0.138196601125011   0.138196601125011; ];
	flag = 0;
end % int_3d_gauss_4

