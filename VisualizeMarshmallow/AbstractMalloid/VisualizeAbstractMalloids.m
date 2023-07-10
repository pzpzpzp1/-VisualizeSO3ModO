clear all; close all;
addpath('Marshmallow')

nonmetric = 0;
dimension = 3;
resolution = 11;
idxs = [1 2 3];

fnameM = ['cache/malloidD' num2str(dimension) 'R' num2str(resolution) '.mat'];
fnameMT = ['cache/malloidTetmeshesD' num2str(dimension) 'R' num2str(resolution) '.mat'];
fnameMAT = ['cache/malloidAbstractTetmeshesD' num2str(dimension) 'R' num2str(resolution) '.mat'];

assert(exist(fnameM)~=0);
assert(numel(idxs)==3);
assert(max(idxs)<= dimension);

load(fnameM);
load(fnameMT);
load(fnameMAT);

%% compute colormaps
Cra = sqrt(sum(malloidAbstractTetmesh.V2P.*malloidAbstractTetmesh.V2P,2));
Caz = atan2(malloidAbstractTetmesh.V2P(:,1),malloidAbstractTetmesh.V2P(:,2));
rescaledCra = (Cra - min(Cra))./max(Cra - min(Cra));
rescaledCaz = (Caz - min(Caz))./max(Caz - min(Caz));

%% base malloid
f = figure; hold on; axis equal; rotate3d on; title('malloid');
scatter3(malloid.V2P(malloid.edgeInds,1),malloid.V2P(malloid.edgeInds,2),malloid.V2P(malloid.edgeInds,3),'b')
scatter3(malloid.V2P(malloid.faceInds,1),malloid.V2P(malloid.faceInds,2),malloid.V2P(malloid.faceInds,3),'g')
scatter3(malloid.V2P(malloid.cornerInds,1),malloid.V2P(malloid.cornerInds,2),malloid.V2P(malloid.cornerInds,3),500,'k.')
PlotMallow(malloid, 0, f, 'green',.1,'-');

%% correctness checks
figure; hold on; axis equal; rotate3d on; title('CORRECTNESS: malloid Tetmesh');
interiorVertInds = 1:size(malloidTetmesh.V2P,1); interiorVertInds([malloidTetmesh.edgeInds malloidTetmesh.faceInds malloidTetmesh.cornerInds])=[];
scatter3(malloidTetmesh.V2P(malloidTetmesh.edgeInds,1),malloidTetmesh.V2P(malloidTetmesh.edgeInds,2),malloidTetmesh.V2P(malloidTetmesh.edgeInds,3),'b')
scatter3(malloidTetmesh.V2P(malloidTetmesh.faceInds,1),malloidTetmesh.V2P(malloidTetmesh.faceInds,2),malloidTetmesh.V2P(malloidTetmesh.faceInds,3),'g')
scatter3(malloidTetmesh.V2P(malloidTetmesh.cornerInds,1),malloidTetmesh.V2P(malloidTetmesh.cornerInds,2),malloidTetmesh.V2P(malloidTetmesh.cornerInds,3),500,'k.')
scatter3(malloidTetmesh.V2P(interiorVertInds,1),malloidTetmesh.V2P(interiorVertInds,2),malloidTetmesh.V2P(interiorVertInds,3),5,'k.')

figure; hold all; axis equal; rotate3d on; title('CORRECTNESS: colorcoded abstract');
scatter3(malloidAbstractTetmesh.V2P(malloidAbstractTetmesh.newFaceInds,1),malloidAbstractTetmesh.V2P(malloidAbstractTetmesh.newFaceInds,2),malloidAbstractTetmesh.V2P(malloidAbstractTetmesh.newFaceInds,3),'g');
scatter3(malloidAbstractTetmesh.V2P(malloidAbstractTetmesh.newEdgeInds,1),malloidAbstractTetmesh.V2P(malloidAbstractTetmesh.newEdgeInds,2),malloidAbstractTetmesh.V2P(malloidAbstractTetmesh.newEdgeInds,3),'b');
scatter3(malloidAbstractTetmesh.V2P(malloidAbstractTetmesh.newCornerInds,1),malloidAbstractTetmesh.V2P(malloidAbstractTetmesh.newCornerInds,2),malloidAbstractTetmesh.V2P(malloidAbstractTetmesh.newCornerInds,3),'k');

%% show base colormap
f = figure; hold all; axis equal; rotate3d on; title('colored radially initial');
scatter3(malloidAbstractTetmesh.V2P(:,1),malloidAbstractTetmesh.V2P(:,2),malloidAbstractTetmesh.V2P(:,3),5,rescaledCra);
PlotMallow(malloid, 0, f, 'green',.05,'none');

f = figure; hold all; axis equal; rotate3d on; title('colored azimuthally initial');
scatter3(malloidAbstractTetmesh.V2P(:,1),malloidAbstractTetmesh.V2P(:,2),malloidAbstractTetmesh.V2P(:,3),5,rescaledCaz);
PlotMallow(malloid, 0, f, 'green',.05,'none');

%% plot MDS embedded
% V = malloidAbstractTetmesh.mds_V2P;
% figure; hold all; axis equal; rotate3d on; title('colored radially embedded mds');
% V = malloidAbstractTetmesh.mds_V2P;
% scatter3(V(:,1),V(:,2),V(:,3),5,rescaledCra);
% 
% figure; hold all; axis equal; rotate3d on; title('colored azimuthally embedded mds');
% V = malloidAbstractTetmesh.mds_V2P;
% scatter3(V(:,1),V(:,2),V(:,3),5,rescaledCaz);
% 
% figure; hold all; axis equal; rotate3d on; title('show edges corners and faces colorcoded embedded mds');
% scatter3(V(:,1),V(:,2),V(:,3),5,'k.');
% scatter3(V(malloidAbstractTetmesh.newFaceInds,1),V(malloidAbstractTetmesh.newFaceInds,2),V(malloidAbstractTetmesh.newFaceInds,3),'g');
% scatter3(V(malloidAbstractTetmesh.newEdgeInds,1),V(malloidAbstractTetmesh.newEdgeInds,2),V(malloidAbstractTetmesh.newEdgeInds,3),'b');
% scatter3(V(malloidAbstractTetmesh.newCornerInds,1),V(malloidAbstractTetmesh.newCornerInds,2),V(malloidAbstractTetmesh.newCornerInds,3),'k');

%% score of cmds
['score is:' num2str(malloidAbstractTetmesh.score)]
V = malloidAbstractTetmesh.cmds_V2P;
EDists = malloidAbstractTetmesh.EDists;
embeddedEdgeDifference = (V(malloidAbstractTetmesh.E2V(:,1),:)-V(malloidAbstractTetmesh.E2V(:,2),:));
embeddedEdgeDistances = sqrt(sum(embeddedEdgeDifference.*embeddedEdgeDifference,2));
idx = sub2ind(size(EDists), malloidAbstractTetmesh.E2V(:,1), malloidAbstractTetmesh.E2V(:,2));
score = norm(EDists(idx)-embeddedEdgeDistances);

%% plot CMDS embedded
V = malloidAbstractTetmesh.cmds_V2P;
figure; hold all; axis equal; rotate3d on; title('colored radially embedded cmds');
scatter3(V(:,idxs(1)),V(:,idxs(2)),V(:,idxs(3)),5,rescaledCra);

figure; hold all; axis equal; rotate3d on; title('colored azimuthally embedded cmds');
scatter3(V(:,idxs(1)),V(:,idxs(2)),V(:,idxs(3)),5,rescaledCaz);

figure; hold all; axis equal; rotate3d on; title('show edges corners and faces colorcoded embedded cmds');
scatter3(V(:,idxs(1)),V(:,idxs(2)),V(:,idxs(3)),5,'k.');
scatter3(V(malloidAbstractTetmesh.newFaceInds,idxs(1)),V(malloidAbstractTetmesh.newFaceInds,idxs(2)),V(malloidAbstractTetmesh.newFaceInds,idxs(3)),'g');
scatter3(V(malloidAbstractTetmesh.newEdgeInds,idxs(1)),V(malloidAbstractTetmesh.newEdgeInds,idxs(2)),V(malloidAbstractTetmesh.newEdgeInds,idxs(3)),'b');
scatter3(V(malloidAbstractTetmesh.newCornerInds,idxs(1)),V(malloidAbstractTetmesh.newCornerInds,idxs(2)),V(malloidAbstractTetmesh.newCornerInds,idxs(3)),500,'k.');

%% animated version
C = Cra;
V0 = malloidAbstractTetmesh.V2P;
V1 = malloidAbstractTetmesh.cmds_V2P;
figure; hold all; axis equal; rotate3d on; title('colored AZ embedded cmds animated');
c = 1;
dt = .04;
while(c < 100); c = c+1;
    for t = [[0:dt:1] fliplr([0:dt:1])]
        cla;
        Vt = V0*(1-t)+V1(:,idxs)*t;
        scatter3(Vt(:,1),Vt(:,2),Vt(:,3),5,C);
        scatter3(Vt(malloidAbstractTetmesh.newEdgeInds,1),Vt(malloidAbstractTetmesh.newEdgeInds,2),Vt(malloidAbstractTetmesh.newEdgeInds,3),'b');
        scatter3(Vt(malloidAbstractTetmesh.newCornerInds,1),Vt(malloidAbstractTetmesh.newCornerInds,2),Vt(malloidAbstractTetmesh.newCornerInds,3),500,'k.');
        
        %scatter3(V(malloidAbstractTetmesh.newFaceInds,idxs(1)),V(malloidAbstractTetmesh.newFaceInds,idxs(2)),V(malloidAbstractTetmesh.newFaceInds,idxs(3)),'g');
        drawnow;
    end
end


%% test score
% dist = graphallshortestpaths(sparse(malloidAbstractTetmesh.EDists));
% V3 = cmdscale(dist,3);
% V6 = cmdscale(dist,6);
% V8 = cmdscale(dist,8);
% V20 = cmdscale(dist,20);
% 
% V = V20;
% EDists = malloidAbstractTetmesh.EDists;
% embeddedEdgeDifference = (V(malloidAbstractTetmesh.E2V(:,1),:)-V(malloidAbstractTetmesh.E2V(:,2),:));
% embeddedEdgeDistances = sqrt(sum(embeddedEdgeDifference.*embeddedEdgeDifference,2));
% idx = sub2ind(size(EDists), malloidAbstractTetmesh.E2V(:,1), malloidAbstractTetmesh.E2V(:,2));
% score = norm(EDists(idx)-embeddedEdgeDistances)
% 
% Vm3 = mdscale(dist,3);
% Vm6 = mdscale(dist,6);
% Vm8 = mdscale(dist,8);
% Vm20 = mdscale(dist,20);
%     

