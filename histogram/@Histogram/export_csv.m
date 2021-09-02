% Tue 24 Nov 16:36:59 +08 2020
function export_csv(obj,filename,c,h)
	fid = fopen(filename,'w');
	if (fid<1)
		error('cannot open file for writing');
	end
	if (nargin()<3)
		c=obj.centre;
	end
	if (nargin<4)
		h = obj.h;
	end
	if (isvector(h)) h = rvec(h); end
	fprintf(fid,'%0.3e;',c); %obj.centre);
	fprintf(fid,'\n');
	for idx=1:size(h,1);
		fprintf(fid,'%0.3e;',h(idx,:));
		fprintf(fid,'\n');
	end
	fclose(fid);
end

