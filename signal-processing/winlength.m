% 2017-02-09 15:52:01.066856764 +0100
%% window length for desired cutoff frequency
%% power at fc is halved
%% H(wf) = 1/sqrt(2) H(f)
%% if the filter window were used as a low pass filter
%% note: the user should prefer a windowed ideal low pass filter
%% TODO, relate this to DOF
function nc = winlength(fc,fs,name)
%	wl.rectwin      = 0.361;
%	wl.triwin       = 0.251;
%	wl.hanwin	= 0.222;
%	wl.kaiserwin	= 0.187;
%	wl.lanczoswin	= 0.098;
	wl.rectwin	= 0.502; 
	wl.triwin	= 0.352; 
	wl.hanwin	= 0.311;
	wl.kaiserwin	= 0.263; 
	wl.lanczoswin	= 0.119;
	nc = (fs/fc)*(1./wl.(name));
end

