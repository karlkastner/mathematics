% Wed Jul 11 18:03:37 MSK 2012
% Karl Kästner, Berlin

% p = 5
function [w b flag] = int_3d_gauss_15()
	w = 6*[
	0.030283678097089
	0.006026785714286
	0.006026785714286
	0.006026785714286
	0.006026785714286
	0.011645249086029
	0.011645249086029
	0.011645249086029
	0.011645249086029
	0.010949141561386
	0.010949141561386
	0.010949141561386
	0.010949141561386
	0.010949141561386
	0.010949141561386];

	b = [
	0.250000000000000   0.250000000000000   0.250000000000000   0.250000000000000
	                0   0.333333333333333   0.333333333333333   0.333333333333333
	0.333333333333333                   0   0.333333333333333   0.333333333333333
	0.333333333333333   0.333333333333333                   0   0.333333333333333
	0.333333333333333   0.333333333333333   0.333333333333333                   0
	0.090909090909091   0.090909090909091   0.090909090909091   0.727272727272727
	0.090909090909091   0.090909090909091   0.727272727272727   0.090909090909091
	0.090909090909091   0.727272727272727   0.090909090909091   0.090909090909091
	0.727272727272727   0.090909090909091   0.090909090909091   0.090909090909091
	0.066550153573664   0.066550153573664   0.433449846426336   0.433449846426336
	0.066550153573664   0.433449846426336   0.066550153573664   0.433449846426336
	0.066550153573664   0.433449846426336   0.433449846426336   0.066550153573664
	0.433449846426336   0.066550153573664   0.066550153573664   0.433449846426336
	0.433449846426336   0.066550153573664   0.433449846426336   0.066550153573664
	0.433449846426336   0.433449846426336   0.066550153573664   0.066550153573664
	];
	flag = 0;
end

