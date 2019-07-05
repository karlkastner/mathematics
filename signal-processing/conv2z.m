% 2015-04-10 11:05:51.081504002 +0200
function A = conv2z(A,f)
	zf = zeros(size(f));
	A = [zf,             zeros(size(f,1),size(A,2)),zf;
             zeros(size(A,1),size(f,2)),A,zeros(size(A,1),size(f,2));
	     zf,zeros(size(f,1),size(A,2)),zf];
	A = conv2(A,f,'same');
	A = A(:,size(f,2)+1:end-size(f,2));
	A = A(size(f,1)+1:end-size(f,1),:);
end

