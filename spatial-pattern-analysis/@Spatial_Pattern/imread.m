% Wed 18 May 13:49:09 CEST 2022
function obj = imread(obj,filename)
	img      = imread(filename);
	if (3 == ndims(img))
		% quick conversion to grayscale
		img = mean(img,3);
	end

	% invert to get biomass proxy
	b = max(img(:))-img;

	n = size(b);

	try
		pgw = csvread([filename(1:end-4),'.pgw']);
		dxy = [pgw(1),abs(pgw(4))];
	catch
		disp('Unable to read pgw file, assuming pixel size of 1 m')
		dxy = 1;
	end

	try

	try
		mskfilename = [filename(1:end-4),'_mask.tif'];
		msk = imread(mskfilename);
	catch
		disp('Unable to load mask, analyzing whole area');
		msk = true(n);
	end

	obj.L   = n.*dxy;
	obj.b   = b;
	obj.msk = msk;
end

