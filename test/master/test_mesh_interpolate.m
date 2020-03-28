% Sat Mar 15 14:14:51 WIB 2014
% Karl Kastner, Berlin

%v.surf.rst allows alos interpolation
% TODO interpolati

% load mesh
%msh = load('msh/test-200.msh');

% qgis :	matplotlib	LinearTriInterpolator

% load depth data
load('mat/dep.mat');
X = dep.X;
Y = dep.Y;
V = dep.depth;

iconst = 0;
%mesh = read_msh('msh/test-200.msh');
mesh = mesh_interpolate(mesh,X,Y,V);
%load('mat/mesh-200.mat','mesh');

px = mesh.P(:,1);
py = mesh.P(:,2);

figure(1);
clf
fdx = find(~isnan(mesh.V));
if (iconst)
	patch(px(mesh.T)',py(mesh.T)',mesh.V','edgecolor',[0.5 0.5 0.5]);
	patch(px(mesh.T(fdx,:))',py(mesh.T(fdx,:))',mesh.V(fdx)','edgecolor','k');
else
	patch(px(mesh.T)',py(mesh.T)',mesh.V(mesh.T)','edgecolor',[0.5 0.5 0.5]);
	patch(px(mesh.T(fdx,:))',py(mesh.T(fdx,:))',mesh.V(mesh.T(fdx,:))','edgecolor','k');
end
caxis([0 50]);
axis equal
hold on;
plot(X,Y,'.k');

figure(2);
clf();
fdx = find(~isnan(mesh.E));
if (iconst)
patch(px(mesh.T)',py(mesh.T)',mesh.E','edgecolor',[0.5 0.5 0.5]);
patch(px(mesh.T(fdx,:))',py(mesh.T(fdx,:))',mesh.E(fdx)','edgecolor','k');
else
	patch(px(mesh.T)',py(mesh.T)',mesh.E(mesh.T)','edgecolor',[0.5 0.5 0.5]);
	patch(px(mesh.T(fdx,:))',py(mesh.T(fdx,:))',mesh.E(mesh.T(fdx,:))','edgecolor','k');
end
axis equal

save('mat/mesh-200.mat','mesh');

if (iconst)
	write_polygon('mesh.shp',px(mesh.T),py(mesh.T),'val',mesh.V);
else
	% fix outliers
	fdx = find(~isnan(mesh.V));
	mesh.V(fdx) = max(0,min(mesh.V(fdx),50));
	v1 = mesh.V(mesh.T(:,1));
	v2 = mesh.V(mesh.T(:,2));
	v3 = mesh.V(mesh.T(:,3));
	vm = 1/3*(v1+v2+v3);
	write_polygon('mesh.shp',px(mesh.T),py(mesh.T),'v1',v1,'v2',v2,'v3',v3,'vm',vm);
end

