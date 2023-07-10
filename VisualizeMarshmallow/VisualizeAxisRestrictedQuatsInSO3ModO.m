clear all; close all;
addpath('Marshmallow');

v = [1,1,1];
n = 100;
malloid = Marshmallow(11,0);

quats = axisRestrictedQuats(v,n);
AA = quatToAxisAngle(quats);
aa = AA(:,1:3).*AA(:,4);

[zs,ds,inds,aaInCenter] = equivalenceInSO3ModO(aa, repmat([0 0 0],size(aa,1),1));
[zs,ds,inds,aaInCenter] = equivalenceInSO3ModO(aaInCenter, repmat([0 0 0],size(aa,1),1));
assert(all(inds))

f = figure; hold on; axis equal; rotate3d on;
PlotMallow(malloid, 0, f, 'green',.1,'none');
patch('Faces',malloid.SE2V,'Vertices',malloid.V2P);
scatter3(aaInCenter(:,1),aaInCenter(:,2),aaInCenter(:,3),'k')







