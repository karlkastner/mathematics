% Wed 18 May 13:49:09 CEST 2022
function [img, alpha, obj] = imread(obj,filename)
	[img, map, alpha] = imread(filename);
	if (~isempty(map))
		error('Convert indexed images to grayscale or RGB first');
	end
	img = double(img);
	if (3 == ndims(img))
		% sloppy conversion to grayscale
		img = mean(img,3);
	end

	% invert to get biomass proxy
	b = max(img(:))-img;
	% stretch contrast to 1
	b = b./max(b(:));

	n = size(b);

	% Try to load pgw from file
	try
		pgw = csvread([filename(1:end-4),'.pgw']);
		obj.stat.pgw = pgw;
		% A    D  B   -E C F
		% wc ws   hs hc
		dxy   = [hypot(pgw(1),pgw(2)),hypot(pgw(3),pgw(4))];
		angle = -pi/2+atan2(pgw(1)/dxy(1),pgw(2)/dxy(2));
		obj.stat.dxy = dxy;
		obj.stat.bangle_rad = angle;
	catch
		disp('Unable to read pgw file, assuming pixel size of 1 m x 1 m')
		dxy = 1;
		obj.stat.pgw = [];
		obj.stat.dxy = NaN;
		obj.stat.bangle_rad = NaN;
	end

	if (~isempty(alpha))
		% note that we can even alias the maks, but we do not do this here
		msk = alpha > 0;
	else
		% try to load mask from file
	try
		mskfilename = [filename(1:end-4),'_mask.tif'];
		msk = imread(mskfilename);
	catch
		disp('Unable to load mask, analyzing whole area');
		msk = true(n);
	end % catch of try
	end % else of ~ isempty alpha

	obj.L   = n.*dxy;
	obj.b   = b;
	obj.msk.b = msk;
end

