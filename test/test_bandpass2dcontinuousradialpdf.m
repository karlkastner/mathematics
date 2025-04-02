L  = 10;
n = 10*L^2;
fr = linspace(0,L,n)';


a = [1/2,1,2];
p = [1/2,1,2];
normalize = true;

for idx=1:length(a)
for jdx=1:length(p)
	S       = bandpass2dcontinuousradialpdf(fr,a(idx),p(jdx),normalize);
	[fc,Sc] = bandpass2dcontinuousradialpdf_mode(a(idx),p(jdx));


subplot(3,3,idx+3*(jdx-1))
cla
plot(fr,S);
hold on
plot(fc,Sc,'*')
xlim([0,3*fc])

II(idx,jdx) = sum(mid(S).*diff(fr));
[aa(idx,jdx),pp(idx,jdx)] = bandpass2dcontinuousradialpdf_mode2par(fc,Sc);



end
end

aa
pp
II


- 2*pi*fc^2/(2*(2*p-1))
