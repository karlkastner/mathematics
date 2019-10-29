% Fri Jul 13 11:59:30 MSK 2012
% Karl Kästner, Berlin

% p = 6
function [w b flag] = int_2d_gauss_24()
	flag = 0;
	w = 6*[
	0.665379170969464506e-2
	0.665379170969464506e-2
	0.665379170969464506e-2
	0.665379170969464506e-2
	0.167953517588677620e-2
	0.167953517588677620e-2
	0.167953517588677620e-2
	0.167953517588677620e-2
	0.9226196923942398432e-2
	0.9226196923942398432e-2
	0.9226196923942398432e-2
	0.9226196923942398432e-2
	0.803571428571428248e-2
	0.803571428571428248e-2
	0.803571428571428248e-2
	0.803571428571428248e-2
	0.803571428571428248e-2
	0.803571428571428248e-2
	0.803571428571428248e-2
	0.803571428571428248e-2
	0.803571428571428248e-2
	0.803571428571428248e-2
	0.803571428571428248e-2
	0.803571428571428248e-2
	];

	b = [
	0.214602871259152   0.214602871259152   0.214602871259152   0.356191386222545
	0.214602871259152   0.214602871259152   0.356191386222545   0.214602871259152
	0.214602871259152   0.356191386222545   0.214602871259152   0.214602871259152
	0.356191386222545   0.214602871259152   0.214602871259152   0.214602871259152
	0.040673958534611   0.040673958534611   0.040673958534611   0.877978124396166
	0.040673958534611   0.040673958534611   0.877978124396166   0.040673958534611
	0.040673958534611   0.877978124396166   0.040673958534611   0.040673958534611
	0.877978124396166   0.040673958534611   0.040673958534611   0.040673958534611
	0.032986329573173   0.322337890142276   0.322337890142276   0.322337890142276
	0.322337890142276   0.032986329573173   0.322337890142276   0.322337890142276
	0.322337890142276   0.322337890142276   0.032986329573173   0.322337890142276
	0.322337890142276   0.322337890142276   0.322337890142276   0.032986329573173
	0.063661001875018   0.063661001875018   0.269672331458316   0.603005664791649
	0.063661001875018   0.063661001875018   0.603005664791649   0.269672331458316
	0.063661001875018   0.269672331458316   0.063661001875018   0.603005664791649
	0.063661001875018   0.269672331458316   0.603005664791649   0.063661001875018
	0.063661001875018   0.603005664791649   0.063661001875018   0.269672331458316
	0.063661001875018   0.603005664791649   0.269672331458316   0.063661001875018
	0.269672331458316   0.063661001875018   0.063661001875018   0.603005664791649
	0.269672331458316   0.063661001875018   0.603005664791649   0.063661001875018
	0.269672331458316   0.603005664791649   0.063661001875018   0.063661001875018
	0.603005664791649   0.063661001875018   0.063661001875018   0.269672331458316
	0.603005664791649   0.063661001875018   0.269672331458316   0.063661001875018
	0.603005664791649   0.269672331458316   0.063661001875018   0.063661001875018
	];
end % int_2d_gauss_24

