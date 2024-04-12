%Tue 20 Feb 14:13:10 CET 2024
function eig_ = eig3(aa)

a = squeeze(aa(1,1,:));
b = squeeze(aa(1,2,:));
c = squeeze(aa(1,3,:));
d = squeeze(aa(2,1,:));
e = squeeze(aa(2,2,:));
f = squeeze(aa(2,3,:));
g = squeeze(aa(3,1,:));
h = squeeze(aa(3,2,:));
i = squeeze(aa(3,3,:));
n = length(a);
o = ones(n,1);
c = [-o, (a + e + i), (b.*d - a.*e - a.*i + c.*g - e.*i + f.*h), a.*e.*i - a.*f.*h - b.*d.*i + b.*f.*g + c.*d.*h - c.*e.*g];
eig_ = roots3(c);
end

