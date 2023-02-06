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
		dxy = [pgw(1),abs(pgw(4))];
	catch
		disp('Unable to read pgw file, assuming pixel size of 1 m x 1 m')
		dxy = 1;
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
	obj.msk = msk;
end

