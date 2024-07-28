% Thu 14 May 19:43:45 +08 2020
% c.f. akin
% Karl Kästner, Berlin
% c.f. Dunavant
% baricentric coordinates and weights for gauss quadrature on the triangle
% 19 points, 9th order accurate
function [w, b, flag] = int_2d_gauss_19()

w =[
	0.097135796282799
	0.031334700227139
	0.031334700227139
	0.031334700227139
	0.077827541004774
	0.077827541004774
	0.077827541004774
	0.079647738927210
	0.079647738927210
	0.079647738927210
	0.025577675658698
	0.025577675658698
	0.025577675658698
	0.043283539377289
	0.043283539377289
	0.043283539377289
	0.043283539377289
	0.043283539377289
	0.043283539377289
];

b = [
0.333333333333333 0.333333333333333 0.333333333333333
0.020634961602525 0.489682519198738 0.489682519198738
0.489682519198738 0.020634961602525 0.489682519198738 
0.489682519198738 0.489682519198738 0.020634961602525
0.125820817014127 0.437089591492937 0.437089591492937
0.437089591492937 0.125820817014127 0.437089591492937
0.437089591492937 0.437089591492937 0.125820817014127
0.623592928761935 0.188203535619033 0.188203535619033
0.188203535619033 0.623592928761935 0.188203535619033
0.188203535619033 0.188203535619033 0.623592928761935
0.910540973211095 0.044729513394453 0.044729513394453
0.044729513394453 0.910540973211095 0.044729513394453
0.044729513394453 0.044729513394453 0.910540973211095
   0.741198598784498   0.221962989160766   0.036838412054736
   0.741198598784498   0.036838412054736   0.221962989160766
   0.221962989160766   0.741198598784498   0.036838412054736
   0.221962989160766   0.036838412054736   0.741198598784498
   0.036838412054736   0.741198598784498   0.221962989160766
   0.036838412054736   0.221962989160766   0.741198598784498
];

end % int_gauss_2d_19
