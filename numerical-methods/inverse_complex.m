% Tue 12 Dec 09:08:34 CET 2017
% Karl Kastner, Berlin

function [re, im] = inverse_complex(re,im)
	a2 = re.^2 + im.^2;
	re = re./a2;
	im = -im./a2;
end

