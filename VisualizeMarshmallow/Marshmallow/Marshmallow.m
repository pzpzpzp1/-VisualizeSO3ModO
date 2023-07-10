% usage: [malloid, malloidTetmesh, malloidAbstractTetmesh]= Marshmallow(5, 1)
function [malloid, malloidTetmesh, malloidAbstractTetmesh]= Marshmallow(resolution, Visualize, dimension)

addpath('AbstractMalloid')

AssertChecks = 0;
digits = 5;
if(nargin ==0)
    resolution = 4;
    Visualize = 0;
    dimension = 3;
end
assert(resolution <= 60)

quality = 1.2;
volumeMax = .0001;


[x,y,z] = sphere;
xrad = resolution; rrad = xrad;
assert(rrad<=xrad);

close all;

% f1 = figure('units','normalized','outerposition',[0 0 1 1]); 
% hold on; axis equal; rotate3d on;
% xlim([-5 5]); ylim([-5 5]); zlim([-5 5])
% ptc = surf(pi*x,pi*y,pi*z,'edgecolor','none','facecolor','green'); alpha(ptc,.3)

R1 = [0 0 1 0]; %Identity
R2 = [0 0 1 pi/2]; %Rx90
R3 = [0 1 0 pi/2]; %Ry90
R4 = [1 1 1 pi/1.5]; % Rxyz120
R5 = [1 0 0 pi/2]; %Rz90

Q1 = axang2quat(R1);
Q2 = axang2quat(R2);
Q3 = axang2quat(R3);
Q4 = axang2quat(R4);
Q5 = axang2quat(R5);

%QS = axang2quat([0 0 1 pi/4]);

plane12 = Q2-Q1; plane12 = plane12 / norm(plane12);
nspX = null(plane12);

plane13 = Q3-Q1; plane13 = plane13 / norm(plane13);
nspY = null(plane13);

plane14 = Q4-Q1; plane14 = plane14 / norm(plane14);
nspXYZ = null(plane14);

plane15 = Q5-Q1; plane15 = plane15 / norm(plane15);
nspZ = null(plane15);

nspXZ = null([plane13; plane12]);
nspXY = null([plane15; plane12]);
nspX_Y_XYZ = null([plane13; plane12; plane14]);
nspX_Z_XYZ = null([plane15; plane12; plane14]);
nspY_Z_XYZ = null([plane15; plane13; plane14]);
nspX_XYZ = null([plane12; plane14]);

AAX = [];
for i = -1:.1:1
    for j = -1:.1:1
        for k = -1:.1:1
            v = nspX*[i;j;k];
            aa = quat2axang(v');
            aa = aa(1:3)*aa(4);
            AAX = [AAX;aa];
        end
    end
end

AAXY = [];
for i = -1:.1:1
    for j = -1:.1:1
            v = nspXY*[i;j];
            aa = quat2axang(v');
            aa = aa(1:3)*aa(4);
            AAXY = [AAXY;aa];
    end
end

AAX_XYZ = [];
for i = -1:.1:1
    for j = -1:.1:1
            v = nspX_XYZ*[i;j];
            aa = quat2axang(v');
            aa = aa(1:3)*aa(4);
            AAX_XYZ = [AAX_XYZ;aa];
    end
end

AAXZ = [];
for i = -1:.1:1
    for j = -1:.1:1
            v = nspXZ*[i;j];
            aa = quat2axang(v');
            aa = aa(1:3)*aa(4);
            AAXZ = [AAXZ;aa];
    end
end

corner1 = quat2axang(nspX_Y_XYZ'); corner1 = corner1(1:3)*corner1(4);
corner2 = quat2axang(nspX_Z_XYZ'); corner2 = corner2(1:3)*corner2(4);
corner3 = quat2axang(nspY_Z_XYZ'); corner3 = corner3(1:3)*corner3(4);

% scatter3(0,0,0,300,'r.')
% scatter3(0,0,pi/2,300,'r.')
% scatter3(0,pi/2,0,300,'r.')
% scatter3(pi*2/3/sqrt(3),pi*2/3/sqrt(3),pi*2/3/sqrt(3),300,'r.')

% scatter3(AAX(:,1),AAX(:,2),AAX(:,3),5,'b')
% scatter3(AAXY(:,1),AAXY(:,2),AAXY(:,3),30,'r')
% scatter3(AAXZ(:,1),AAXZ(:,2),AAXZ(:,3),30,'r')
% scatter3(AAX_XYZ(:,1),AAX_XYZ(:,2),AAX_XYZ(:,3),30,'r')
% scatter3(corner1(1),corner1(2),corner1(3),1000,'k.')
% scatter3(corner2(1),corner2(2),corner2(3),1000,'k.')

%% confirmed corner1 and corner2 are equal to eachothother but rotated a third about 1,1,1.

res = 100;
c1toc2 = [linspace(nspX_Y_XYZ(1,:),nspX_Z_XYZ(1,:),res);...
    linspace(nspX_Y_XYZ(2,:),nspX_Z_XYZ(2,:),res);...
    linspace(nspX_Y_XYZ(3,:),nspX_Z_XYZ(3,:),res);...
    linspace(nspX_Y_XYZ(4,:),nspX_Z_XYZ(4,:),res)];
c1toc2AA = quat2axang(c1toc2'); c1toc2AA = c1toc2AA(:,1:3).*c1toc2AA(:,4);

c1toc2AAe1 = c1toc2AA;
c1toc2AAe2 = c1toc2AA*axang2rotm([1 1 1 pi/1.5])';
c1toc2AAe3 = c1toc2AA*axang2rotm([1 1 1 pi*2/1.5])';

triangleSurfaceX = []; resX = xrad;
B1 = linspace(0,1,resX);
for b1i = 1:numel(B1)
    b1 = B1(b1i);
    resYi = resolution - b1i + 1; %ceil((numel(B1)-b1i+1)/numel(B1)*resX);
    for b2= linspace(0,(1-b1),resYi)
        b3= 1-b1-b2;
        v = [nspY_Z_XYZ nspX_Z_XYZ nspX_Y_XYZ]*[b1;b2;b3];
        aa = quat2axang(v');
        aa = aa(1:3)*aa(4);
        triangleSurfaceX = [triangleSurfaceX; aa];
    end
end

faceCenterQ = axang2quat([1 0 0 pi/4])';
cornerwedgeSurfaceX = []; resR = resX; resX = rrad; 
B1 = linspace(0,1,resX);
for b1i = 1:numel(B1)
    b1 = B1(b1i);
    resYi = resolution - b1i + 1; %ceil((numel(B1)-b1i+1)/numel(B1)*resR);
    for b2= linspace(0,(1-b1),resYi)
        b3= 1-b1-b2;
        v = [faceCenterQ nspY_Z_XYZ nspX_Z_XYZ ]*[b1;b2;b3];
        aa = quat2axang(v');
        aa = aa(1:3)*aa(4);
        cornerwedgeSurfaceX = [cornerwedgeSurfaceX; aa];
    end
end

aa =  quat2axang(nspX_Z_XYZ'); aa = aa(1:3)*aa(4); aa = axang2rotm([1 0 0 -pi/2])*aa'; otherCorner = axang2quat([aa' norm(aa)])';
widewedgeSurfaceX = []; resX = resX; resR = resR;
B1 = linspace(0,1,resX);
for b1i = 1:numel(B1)
    b1 = B1(b1i);
    resYi = resolution - b1i + 1; %ceil((numel(B1)-b1i+1)/numel(B1)*resR);
    for b2= linspace(0,(1-b1),resYi)
        b3= 1-b1-b2;
        v = [faceCenterQ nspY_Z_XYZ otherCorner]*[b1;b2;b3];
        aa = quat2axang(v');
        aa = aa(1:3)*aa(4);
        widewedgeSurfaceX = [widewedgeSurfaceX; aa];
    end
end

% f2 = figure('units','normalized','outerposition',[0 0 1 1]); 
% hold on; axis equal; rotate3d on;
%xlim([-5 5]); ylim([-5 5]); zlim([-5 5])
%ptc = surf(pi*x,pi*y,pi*z,'edgecolor','none','facecolor','green'); alpha(ptc,.3)
% 
% plot3(c1toc2AAe1(:,1),c1toc2AAe1(:,2),c1toc2AAe1(:,3),'r','LineWidth',2)
% plot3(c1toc2AAe2(:,1),c1toc2AAe2(:,2),c1toc2AAe2(:,3),'r','LineWidth',2)
% % plot3(c1toc2AAe3(:,1),c1toc2AAe3(:,2),c1toc2AAe3(:,3),'r','LineWidth',2)

if(Visualize)
    figure; hold on; rotate3d on; axis equal;
    scatter3(triangleSurfaceX(:,1),triangleSurfaceX(:,2),triangleSurfaceX(:,3),'r');
    scatter3(cornerwedgeSurfaceX(:,1),cornerwedgeSurfaceX(:,2),cornerwedgeSurfaceX(:,3),'k');
    scatter3(widewedgeSurfaceX(:,1),widewedgeSurfaceX(:,2),widewedgeSurfaceX(:,3),'b');
end

% % fun fact, widewedgeSurfaceX and cornerwedgeSurfaceX are the same but
% % rotated by 2pi/8.
% s2 = (axang2rotm([1,0,0,-pi/4])*cornerwedgeSurfaceX')';
% scatter3(s2(:,1),s2(:,2),s2(:,3),20,'k');
% combinedV = [s2; widewedgeSurfaceX];
% [uniqueV, ia, ib] = unique(fix(combinedV*(10^digits)),'rows');
% repeatInds = find(accumarray(ib,1)>1);

% scatter3(0,0,0,300,'r.')
% scatter3(0,0,pi/2,300,'r.')
% scatter3(0,pi/2,0,300,'r.')
% scatter3(pi*2/3/sqrt(3),pi*2/3/sqrt(3),pi*2/3/sqrt(3),300,'r.')
% ptc = surf(pi*x/4,pi*y/4,pi*z/4,'edgecolor','none','facecolor','green'); alpha(ptc,.3)
% 
% plot3(triangleSurfaceX(:,1),triangleSurfaceX(:,2),triangleSurfaceX(:,3),'r');
% plot3(cornerwedgeSurfaceX(:,1),cornerwedgeSurfaceX(:,2),cornerwedgeSurfaceX(:,3),'k');
% plot3(widewedgeSurfaceX(:,1),widewedgeSurfaceX(:,2),widewedgeSurfaceX(:,3),'b');



%% surface vertices created. Time to make the triangle mesh.
% parse triangleSurfaceX as xrad X xrad
% parse cornerwedgeSurfaceX as xrad X rrad
% parse widewedgeSurfaceX as xrad X rrad

T = [];
Vinds = 1:size(triangleSurfaceX,1);
assert(sum(1:resolution) == numel(Vinds));
rownum = 1;
while numel(Vinds)~=0
    row1 = Vinds(1:xrad-rownum+1);
    row2 = Vinds(xrad-rownum+2:2*xrad-2*rownum+1);
    assert(numel(row2)+1==numel(row1));
    
    for i = 1:numel(row2)
        T = [T;row1(i) row1(i+1) row2(i)];
        if(i~=numel(row2))
            T = [T;row2(i) row1(i+1) row2(i+1)];
        end
    end
    
    Vinds = Vinds(xrad-rownum+2:end);
    rownum = rownum+1;
end
T1 = T;

T = [];
Vinds = 1:size(cornerwedgeSurfaceX,1);
assert(sum(1:resolution) == numel(Vinds));
rownum = 1;
while numel(Vinds)~=0
    row1 = Vinds(1:xrad-rownum+1);
    row2 = Vinds(xrad-rownum+2:2*xrad-2*rownum+1);
    assert(numel(row2)+1==numel(row1));
    
    for i = 1:numel(row2)
        T = [T;row1(i) row1(i+1) row2(i)];
        if(i~=numel(row2))
            T = [T;row2(i) row1(i+1) row2(i+1)];
        end
    end
    
    Vinds = Vinds(xrad-rownum+2:end);
    rownum = rownum+1;
end
T2 = T;

T = [];
Vinds = 1:size(widewedgeSurfaceX,1);
assert(sum(1:resolution) == numel(Vinds));
rownum = 1;
while numel(Vinds)~=0
    row1 = Vinds(1:xrad-rownum+1);
    row2 = Vinds(xrad-rownum+2:2*xrad-2*rownum+1);
    assert(numel(row2)+1==numel(row1));
    
    for i = 1:numel(row2)
        T = [T;row1(i) row1(i+1) row2(i)];
        if(i~=numel(row2))
            T = [T;row2(i) row1(i+1) row2(i+1)];
        end
    end
    
    Vinds = Vinds(xrad-rownum+2:end);
    rownum = rownum+1;
end
T3 = T;


mesh1.V2P = triangleSurfaceX;
mesh2.V2P = cornerwedgeSurfaceX;
mesh3.V2P = widewedgeSurfaceX;
mesh1.T2V = T1;
mesh2.T2V = T2;
mesh3.T2V = T3;

mc1 = mesh1;
mc2 = mesh1; mc2.V2P = (axang2rotm([1 0 0 pi/2])*mc1.V2P')';
mc3 = mesh1; mc3.V2P = (axang2rotm([1 0 0 pi])*mc1.V2P')';
mc4 = mesh1; mc4.V2P = (axang2rotm([1 0 0 -pi/2])*mc1.V2P')';
mc5 = mesh1; mc5.V2P = (axang2rotm([0 0 1 pi/2])*mc1.V2P')';
mc6 = mc5; mc6.V2P = (axang2rotm([1 0 0 pi/2])*mc5.V2P')';
mc7 = mc5; mc7.V2P = (axang2rotm([1 0 0 pi])*mc5.V2P')';
mc8 = mc5; mc8.V2P = (axang2rotm([1 0 0 -pi/2])*mc5.V2P')';

mf = CombineMesh(mesh2, mesh3,1);
mq2 = mf; mq2.V2P = (axang2rotm([1 0 0 pi/2])*mq2.V2P')';
mf = CombineMesh(mf, mq2,1);
mq2 = mf; mq2.V2P = (axang2rotm([1 0 0 pi])*mq2.V2P')';
mf = CombineMesh(mf, mq2,1);
facemesh = mf;

%% mf is now a full face.
mq2 = mf; mq2.V2P = (axang2rotm([0 1 0 pi/2])*mq2.V2P')';
mf = CombineMesh(mf, mq2,1);
mq2 = mf; mq2.V2P = (axang2rotm([0 1 0 pi])*mq2.V2P')';
mf = CombineMesh(mf, mq2,1);

mq2 = facemesh; mq2.V2P = (axang2rotm([0 0 1 pi/2])*mq2.V2P')';
mq3 = facemesh; mq3.V2P = (axang2rotm([0 0 1 -pi/2])*mq3.V2P')';
mf = CombineMesh(mf, mq2,1);
mf = CombineMesh(mf, mq3,1);

%% add corners now
mf = CombineMesh(mf, mc1,1);
mf = CombineMesh(mf, mc2,1);
mf = CombineMesh(mf, mc3,1);
mf = CombineMesh(mf, mc4,1);
mf = CombineMesh(mf, mc5,1);
mf = CombineMesh(mf, mc6,1);
mf = CombineMesh(mf, mc7,1);
mf = CombineMesh(mf, mc8,1);

malloid = mf;

if(Visualize)
    figure; hold on; rotate3d on; axis equal;
    ptc = patch('Faces',malloid.T2V,'Vertices',malloid.V2P,'FaceColor','green','LineStyle','-'); alpha(ptc,.1);
end

[V1,Qout] = shiftPoints(malloid.V2P, axang2quat([0,0,1,pi/2]));
[V2,Qout] = shiftPoints(malloid.V2P, axang2quat([0,1,0,pi/2]));
[V3,Qout] = shiftPoints(malloid.V2P, axang2quat([1,0,0,pi/2]));
[V4,Qout] = shiftPoints(malloid.V2P, axang2quat([0,0,-1,pi/2]));
[V5,Qout] = shiftPoints(malloid.V2P, axang2quat([0,-1,0,pi/2]));
[V6,Qout] = shiftPoints(malloid.V2P, axang2quat([-1,0,0,pi/2]));

[VC1,Qout] = shiftPoints(malloid.V2P, axang2quat([1,1,1,2*pi/3]));
[VC2,Qout] = shiftPoints(malloid.V2P, axang2quat([1,1,-1,2*pi/3]));
[VC3,Qout] = shiftPoints(malloid.V2P, axang2quat([1,-1,1,2*pi/3]));
[VC4,Qout] = shiftPoints(malloid.V2P, axang2quat([1,-1,-1,2*pi/3]));
[VC5,Qout] = shiftPoints(malloid.V2P, axang2quat([-1,1,1,2*pi/3]));
[VC6,Qout] = shiftPoints(malloid.V2P, axang2quat([-1,1,-1,2*pi/3]));
[VC7,Qout] = shiftPoints(malloid.V2P, axang2quat([-1,-1,1,2*pi/3]));
[VC8,Qout] = shiftPoints(malloid.V2P, axang2quat([-1,-1,-1,2*pi/3]));

combinedV = [malloid.V2P;V1;V2;V3;V4;V5;V6;VC1;VC2;VC3;VC4;VC5;VC6;VC7;VC8];
%% properties of combinedV. original faces are duplicated twice. original edges are duplicated 4 times.

[uniqueV, ia, ib] = unique(fix(combinedV*(10^digits)),'rows');
nV = size(malloid.V2P,1);
onlyOnceInds = find(accumarray(ib,1)==1);
onlyTwiceInds = find(accumarray(ib,1)==2);
onlyThreeXInds = find(accumarray(ib,1)==3);
onlyFourXInds = find(accumarray(ib,1)==4);

if(0)
    figure; hold on; rotate3d on; title('make sure the right verts are matched');
    scatter3(uniqueV(:,1),uniqueV(:,2),uniqueV(:,3))
    scatter3(uniqueV(onlyOnceInds,1),uniqueV(onlyOnceInds,2),uniqueV(onlyOnceInds,3),'r')
    scatter3(uniqueV(onlyTwiceInds,1),uniqueV(onlyTwiceInds,2),uniqueV(onlyTwiceInds,3),'g')
    scatter3(uniqueV(onlyThreeXInds,1),uniqueV(onlyThreeXInds,2),uniqueV(onlyThreeXInds,3),'b')
    scatter3(uniqueV(onlyFourXInds,1),uniqueV(onlyFourXInds,2),uniqueV(onlyFourXInds,3),'k.')
end

repeatInds = find(accumarray(ib,1)==4);
allinds = 1:size(combinedV,1);
indsOfUniqueV = allinds(ia);
edgeInds = indsOfUniqueV(onlyThreeXInds);
edgeInds(edgeInds > nV)=[];
cornerInds = indsOfUniqueV(onlyFourXInds);
cornerInds(cornerInds> nV)=[];
faceInds = indsOfUniqueV(onlyTwiceInds);
faceInds(faceInds> nV)=[];

% ASSIGN THE IMPORTANT FEATURES!
malloid.faceInds = faceInds;
malloid.cornerInds = cornerInds;
malloid.edgeInds = edgeInds;

% remains to identify vertices to eachother.
%ib = ib(1:nV);
inds = repmat(1:nV,1,15);
spr = sparse(ib,1:size(ib,1),1:size(ib,1));
faceindsb = find(sum(spr>0,2)==2);
sprf = 0*spr; sprf(faceindsb,:)=spr(faceindsb,:);
assert(all(sum(sprf>0,2)==2 | sum(sprf>0,2)==0))
[ii,jj] = find(sprf);
[~, sortind] = sort(ii);
jj = jj(sortind);
pairinds = reshape(jj,2,[])';
pairinds = unique(sort(inds(pairinds),2),'rows');

edgeindsb = find(sum(spr>0,2)==3);
sprf = 0*spr; sprf(edgeindsb,:)=spr(edgeindsb,:);
[ii,jj] = find(sprf);
[~, sortind] = sort(ii);
jj = jj(sortind);
tripleinds = reshape(jj,3,[])';
tripleinds = unique(sort(inds(tripleinds),2),'rows');

cornerindsb = find(sum(spr>0,2)==4);
sprf = 0*spr; sprf(cornerindsb,:)=spr(cornerindsb,:);
[ii,jj] = find(sprf);
[~, sortind] = sort(ii);
jj = jj(sortind);
quadinds = reshape(jj,4,[])';
quadinds = unique(sort(inds(quadinds),2),'rows');

% I don't understand why but pair inds contains triple inds, and triple
% inds contains quadinds. gotta prune them.

prune = ismember(pairinds, faceInds(:));
assert(all(sum(prune,2)==0 | sum(prune,2)==2));
pairinds(find(sum(prune,2)==0),:)=[];

prune = ismember(tripleinds, edgeInds);
assert(all(sum(prune,2)==0 | sum(prune,2)==3));
tripleinds(sum(prune,2)==0,:)=[];

prune = ismember(quadinds, cornerInds);
assert(all(sum(prune,2)==0 | sum(prune,2)==4));
quadinds(sum(prune,2)==0,:)=[];

% verify correctness
assert(all(sort(pairinds(:))==sort(faceInds)'))
assert(all(sort(tripleinds(:))==sort(edgeInds)'))
assert(all(sort(quadinds(:))==sort(cornerInds)'))
assert(size(quadinds,1)==6); % there's only 6 unique corners in this shape regardless of resolution.

% verify identifications
if(AssertChecks)
    for i = 1:size(pairinds,1)
        pair = pairinds(i,:);
        assert(equivalenceInSO3ModO(malloid.V2P(pair(1),:), malloid.V2P(pair(2),:)));
    end

    for i = 1:size(tripleinds,1)
        triple = tripleinds(i,:);
        assert(equivalenceInSO3ModO(malloid.V2P(triple(1),:), malloid.V2P(triple(2),:)));
        assert(equivalenceInSO3ModO(malloid.V2P(triple(3),:), malloid.V2P(triple(2),:)));
    end

    for i = 1:size(quadinds,1)
        fourInds = quadinds(i,:);
        assert(equivalenceInSO3ModO(malloid.V2P(fourInds(1),:), malloid.V2P(fourInds(2),:)));
        assert(equivalenceInSO3ModO(malloid.V2P(fourInds(3),:), malloid.V2P(fourInds(2),:)));
        assert(equivalenceInSO3ModO(malloid.V2P(fourInds(3),:), malloid.V2P(fourInds(4),:)));
    end
end

if(Visualize)
    figure; hold on; axis equal; rotate3d on; title('faces edges corners of surface mesh');
    scatter3(malloid.V2P(edgeInds,1),malloid.V2P(edgeInds,2),malloid.V2P(edgeInds,3),'b')
    scatter3(malloid.V2P(faceInds,1),malloid.V2P(faceInds,2),malloid.V2P(faceInds,3),'g')
    scatter3(malloid.V2P(cornerInds,1),malloid.V2P(cornerInds,2),malloid.V2P(cornerInds,3),500,'k.')
end

assert(all(sort([faceInds, edgeInds, cornerInds])==1:size(malloid.V2P,1)));

if nargout > 1
    fid = fopen('Convert/malloid.smesh','w');
    SaveMalloidAsSmesh(fid,malloid);
    fclose(fid);
    
    % add additional points
    %     fid = fopen('Convert/malloid-a.node','w');
    %     nInsert = (3*resolution)^3;
    %     fprintf(fid,'%d 3 0 0\n',nInsert);
    %     Ps = [[1:nInsert]; (rand(3,nInsert)-.5)*pi/2*.99];
    %     fprintf(fid,'%d %f %f %f\n',Ps(:));
    %     fclose(fid);
    
    % create tet mesh. 
    % [a,b,c] = Marshmallow(3,4);
    if(ispc)
        [status,cmdout] = system(['.\..\..\tetgen1.5.1-beta1\build\tetgen.exe .\Convert\malloid.smesh -Ypq' num2str(quality) 'a' num2str(volumeMax)]);
    else
        [status,cmdout] = system(['./../../tetgen1.5.1-beta1/build/tetgen ./Convert/malloid.smesh -Ypq' num2str(quality) 'a' num2str(volumeMax)]);
    end
    assert(status == 0);
    
    % read tet mesh.
    fid = fopen('Convert/malloid.1.ele','r');
    nTets = fscanf(fid,'%d 4 0/n',1);
    tets = fscanf(fid,'%d',5*nTets); tets = reshape(tets,5,nTets)'; tets(:,1) = [];
    fclose(fid);
    
    fid = fopen('Convert/malloid.1.node','r');
    nVerts = fscanf(fid,'%d 3 0 0/n',1);
    verts = fscanf(fid,'%f',4*nVerts); verts = reshape(verts,4,nVerts)'; verts(:,1) = [];
    fclose(fid);
    
    malloidTetmesh.T2V = tets;
    malloidTetmesh.V2P = verts;
    malloidTetmesh.E2V = unique(sort(reshape(malloidTetmesh.T2V(:,[1,2,2,3,3,4,1,3,1,4,2,4])',2,[])',2),'rows');
    
    if(0)
        hold on;
        %scatter3(malloidTetmesh.V2P(:,1),malloidTetmesh.V2P(:,2),malloidTetmesh.V2P(:,3),5,'k')
        %scatter3(Ps(2,:),Ps(3,:),Ps(4,:),5,'k')
        for i = 1:size(malloidTetmesh.E2V, 1)
            vs = malloidTetmesh.V2P(malloidTetmesh.E2V(i,:),:);
            plot3(vs(:,1),vs(:,2),vs(:,3),'k');
        end
    end
    
    % this makes sure tet mesh faces corners and edges still have same
    % indices
    assert(norm(malloidTetmesh.V2P([edgeInds faceInds cornerInds],:)-malloid.V2P([edgeInds faceInds cornerInds],:))<.00001);
    malloidTetmesh.edgeInds = malloid.edgeInds;
    malloidTetmesh.faceInds = malloid.faceInds;
    malloidTetmesh.cornerInds = malloid.cornerInds;
    
    if(Visualize)
        figure; hold on; axis equal; rotate3d on; title('face edges and corners of Tetmesh');
        scatter3(malloidTetmesh.V2P(edgeInds,1),malloidTetmesh.V2P(edgeInds,2),malloidTetmesh.V2P(edgeInds,3),'b')
        scatter3(malloidTetmesh.V2P(faceInds,1),malloidTetmesh.V2P(faceInds,2),malloidTetmesh.V2P(faceInds,3),'g')
        scatter3(malloidTetmesh.V2P(cornerInds,1),malloidTetmesh.V2P(cornerInds,2),malloidTetmesh.V2P(cornerInds,3),500,'k.')
    end
    
end

malloid.pairInds = pairinds;
malloid.tripleInds = tripleinds;
malloid.quadInds = quadinds;

%% for malloidTetmesh, create tetmesh, match vertices, compute edge distances.
if nargout > 2
    newFaceInds = malloid.faceInds;
    newEdgeInds = malloid.edgeInds;
    newCornerInds = malloid.cornerInds;
    
    malloidAbstractTetmesh = malloidTetmesh;
    malloidAbstractTetmesh = rmfield(malloidAbstractTetmesh, 'faceInds');
    malloidAbstractTetmesh = rmfield(malloidAbstractTetmesh, 'edgeInds');
    malloidAbstractTetmesh = rmfield(malloidAbstractTetmesh, 'cornerInds');
    
    % identify vertices based on malloid indices.
    % identify faces
    matchedInds = pairinds;
    Vs = 1:size(malloidAbstractTetmesh.V2P,1);
    Vs(matchedInds(:,2)) = matchedInds(:,1);
    malloidAbstractTetmesh.T2V = Vs(malloidAbstractTetmesh.T2V);

    verts2remove = matchedInds(:,2);
    Vs = 1:size(malloidAbstractTetmesh.V2P,1);
    Vs(verts2remove)=[];
    V2s = 1:numel(Vs);
    V3s = 1:numel(Vs);
    V3s(Vs)=V2s;
    malloidAbstractTetmesh.V2P(verts2remove,:)=[];
    malloidAbstractTetmesh.T2V = V3s(malloidAbstractTetmesh.T2V);
    quadinds = V3s(quadinds);
    tripleinds = V3s(tripleinds);
    newFaceInds(ismember(newFaceInds,verts2remove))=[];
    newEdgeInds(ismember(newEdgeInds,verts2remove))=[];
    newCornerInds(ismember(newCornerInds,verts2remove))=[];
    newFaceInds = V3s(newFaceInds);
    newEdgeInds = V3s(newEdgeInds);
    newCornerInds = V3s(newCornerInds);
    
    % identify corners
    matchedInds = tripleinds;
    Vs = 1:size(malloidAbstractTetmesh.V2P,1);
    Vs(matchedInds(:,2)) = matchedInds(:,1);
    Vs(matchedInds(:,3)) = matchedInds(:,1);
    malloidAbstractTetmesh.T2V = Vs(malloidAbstractTetmesh.T2V);

    verts2remove = matchedInds(:,2:3);
    Vs = 1:size(malloidAbstractTetmesh.V2P,1);
    Vs(verts2remove)=[];
    V2s = 1:numel(Vs);
    V3s = 1:numel(Vs);
    V3s(Vs)=V2s;
    malloidAbstractTetmesh.V2P(verts2remove,:)=[];
    malloidAbstractTetmesh.T2V = V3s(malloidAbstractTetmesh.T2V);
    quadinds = V3s(quadinds);
    newFaceInds(ismember(newFaceInds,verts2remove))=[];
    newEdgeInds(ismember(newEdgeInds,verts2remove))=[];
    newCornerInds(ismember(newCornerInds,verts2remove))=[];
    newFaceInds = V3s(newFaceInds);
    newEdgeInds = V3s(newEdgeInds);
    newCornerInds = V3s(newCornerInds);
    
    % identify corners
    matchedInds = quadinds;
    Vs = 1:size(malloidAbstractTetmesh.V2P,1);
    Vs(matchedInds(:,2)) = matchedInds(:,1);
    Vs(matchedInds(:,3)) = matchedInds(:,1);
    Vs(matchedInds(:,4)) = matchedInds(:,1);
    malloidAbstractTetmesh.T2V = Vs(malloidAbstractTetmesh.T2V);

    verts2remove = matchedInds(:,2:4);
    Vs = 1:size(malloidAbstractTetmesh.V2P,1);
    Vs(verts2remove)=[];
    V2s = 1:numel(Vs);
    V3s = 1:numel(Vs);
    V3s(Vs)=V2s;
    malloidAbstractTetmesh.V2P(verts2remove,:)=[];
    malloidAbstractTetmesh.T2V = V3s(malloidAbstractTetmesh.T2V);
    newFaceInds(ismember(newFaceInds,verts2remove))=[];
    newEdgeInds(ismember(newEdgeInds,verts2remove))=[];
    newCornerInds(ismember(newCornerInds,verts2remove))=[];
    newFaceInds = V3s(newFaceInds);
    newEdgeInds = V3s(newEdgeInds);
    newCornerInds = V3s(newCornerInds);
    
    % reassign the new faceinds edgeinds cornerinds
    malloidAbstractTetmesh.newFaceInds = newFaceInds;
    malloidAbstractTetmesh.newEdgeInds = newEdgeInds;
    malloidAbstractTetmesh.newCornerInds = newCornerInds;
    
    assert(numel(malloidAbstractTetmesh.newCornerInds)==6);
    assert(numel(malloidAbstractTetmesh.newEdgeInds)/(resolution-2)==12);
    assert(numel(malloidAbstractTetmesh.newEdgeInds)/(resolution-2)==12);
    nVinWedge = (resolution-1)*resolution/2;
    nVertsInOctFaceAlone = nVinWedge*8 - 8*(resolution-1) + 1;
    nVertsInTriFaceAlone = (resolution-3)*(resolution-2)/2;
    assert(numel(malloidAbstractTetmesh.newFaceInds) == (nVertsInOctFaceAlone*6+nVertsInTriFaceAlone*8)/2);
    
    if(Visualize)
        figure; hold all; axis equal; rotate3d on; title('show edges corners and faces colorcoded initial');
        V = malloidAbstractTetmesh.V2P;
        %plot3(V(:,1),V(:,2),V(:,3),'.');
        scatter3(V(newFaceInds,1),V(newFaceInds,2),V(newFaceInds,3),'g');
        scatter3(V(newEdgeInds,1),V(newEdgeInds,2),V(newEdgeInds,3),'b');
        scatter3(V(newCornerInds,1),V(newCornerInds,2),V(newCornerInds,3),'k');
        clear('V');
    end
    
    % create edges E2V (sorted of course)
    malloidAbstractTetmesh.E2V = unique(sort(reshape(malloidAbstractTetmesh.T2V(:,[1,2,2,3,3,4,1,3,1,4,2,4])',2,[])',2),'rows');
    
    % create sparse distances
    n = size(malloidAbstractTetmesh.V2P,1);
    'calculating edists'
    tic
    EDists = sparse(n,n);
    [~,ds] = equivalenceInSO3ModO(malloidAbstractTetmesh.V2P(malloidAbstractTetmesh.E2V(:,1),:), malloidAbstractTetmesh.V2P(malloidAbstractTetmesh.E2V(:,2),:));
    EDists(sub2ind(size(EDists), malloidAbstractTetmesh.E2V(:,1), malloidAbstractTetmesh.E2V(:,2))) = ds;
    EDists = EDists + EDists';
    malloidAbstractTetmesh.EDists = EDists;
    toc
%     hist(EDists(find(EDists~=0)));
%     pause
    
    % compute mds embedding
    'calculating all shortest paths'
    tic
    dist = graphallshortestpaths(sparse(EDists));
    toc
    
    %malloidAbstractTetmesh.AllDists = dist;

     %% FULL GLOBAL DISTANCES: TOO SLOW TO USE EVEN WHEN VECTORIZED.    
%     nV1 = size(malloidAbstractTetmesh.V2P,1);
%     i1 = repmat(1:nV1,1,nV1)';
%     i2 = repmat(1:nV1,nV1,1); i2 = i2(:);
%     [~,ds] = equivalenceInSO3ModO(malloidAbstractTetmesh.V2P(i1,:),malloidAbstractTetmesh.V2P(i2,:));
%     dist = zeros(nV1);
%     dist(sub2ind(size(EDists), i1, i2)) = ds;
    
%     'calculating mdscale'
%     tic
%     malloidAbstractTetmesh.mds_V2P = mdscale(dist,dimension);
%     toc
    'calculating cmdscale'
    tic
    malloidAbstractTetmesh.cmds_V2P = cmdscale(dist,dimension);
    toc
    'done computing embedding'
    
%     'calculating mdscale on just EDistsLocal'
%     tic
%     EDistsFull = full(EDists); EDistsFull(find(EDistsFull == 0)) = 0/0;
%     EDistsFull(1:size(EDists,2)+1:end)=0;
%     malloidAbstractTetmesh.mdsLocal_V2P = mdscale(EDistsFull,dimension,'Start','random');
%     toc
%     'done computing embeddings'
    
    %% check how accurate distances are in the embedding cmds
    'calculating score'
    tic
    V = malloidAbstractTetmesh.cmds_V2P;
    embeddedEdgeDifference = (V(malloidAbstractTetmesh.E2V(:,1),:)-V(malloidAbstractTetmesh.E2V(:,2),:));
    embeddedEdgeDistances = sqrt(sum(embeddedEdgeDifference.*embeddedEdgeDifference,2));
    idx = sub2ind(size(EDists), malloidAbstractTetmesh.E2V(:,1), malloidAbstractTetmesh.E2V(:,2));
    score = norm(EDists(idx)-embeddedEdgeDistances);
    malloidAbstractTetmesh.score = score;
    toc
    
    if(Visualize)
        figure; hold all; axis equal; rotate3d on; title('show edges corners and faces colorcoded embedded');
        scatter3(V(malloidAbstractTetmesh.newFaceInds,1),V(malloidAbstractTetmesh.newFaceInds,2),V(malloidAbstractTetmesh.newFaceInds,3),'g');
        scatter3(V(newEdgeInds,1),V(newEdgeInds,2),V(newEdgeInds,3),'b');
        scatter3(V(newCornerInds,1),V(newCornerInds,2),V(newCornerInds,3),'k');
    end
    
    % compute colormap
    Cra = sqrt(sum(malloidAbstractTetmesh.V2P.*malloidAbstractTetmesh.V2P,2));
    Caz = atan2(malloidAbstractTetmesh.V2P(:,1),malloidAbstractTetmesh.V2P(:,2));
    rescaledCra = (Cra - min(Cra))./max(Cra - min(Cra));
    rescaledCaz = (Caz - min(Caz))./max(Caz - min(Caz));
    
    if(Visualize)
        % visualize colormaps on initial and embedded spaces
        f = figure; hold all; axis equal; rotate3d on; title('colored radially initial');
        scatter3(malloidAbstractTetmesh.V2P(:,1),malloidAbstractTetmesh.V2P(:,2),malloidAbstractTetmesh.V2P(:,3),5,rescaledCra);
        PlotMallow(malloid, 0, f, 'green',.1,'none');
        
        f = figure; hold all; axis equal; rotate3d on; title('colored azimuthally initial');
        scatter3(malloidAbstractTetmesh.V2P(:,1),malloidAbstractTetmesh.V2P(:,2),malloidAbstractTetmesh.V2P(:,3),5,rescaledCaz);
        PlotMallow(malloid, 0, f, 'green',.1,'none');
        
        figure; hold all; axis equal; rotate3d on; title('colored radially embedded');
        scatter3(V(:,1),V(:,2),V(:,3),5,rescaledCra);
        figure; hold all; axis equal; rotate3d on; title('colored azimuthally embedded');
        scatter3(V(:,1),V(:,2),V(:,3),5,rescaledCaz);
    end
    
%     for i=1:size(malloidAbstractTetmesh.E2V,1)
%         ei = malloidAbstractTetmesh.E2V(i,:);
%         vs = V(ei,:);
%         plot3(vs(:,1),vs(:,2),vs(:,3));
%     end

end

%% add malloid edge edges alone
E2V = unique(reshape(malloid.T2V(:,[1 2 1 3 2 3])',2,[])','rows');
malloid.E2V = E2V;
keepEdges = ismember(E2V, [malloid.cornerInds malloid.edgeInds]);
keepEdges = keepEdges(:,1) & keepEdges(:,2);
E2V(find(~keepEdges),:)=[];
malloid.SE2V = E2V; % sharp edges. the ones between cornerInds and edgeInds

rightOctagonVertInds = find(malloid.V2P(:,1) > .9*pi/4); 
all(ismember(malloid.T2V,rightOctagonVertInds),2);
rightOctagonTriInds = find(all(ismember(malloid.T2V,rightOctagonVertInds),2));
malloid.rightOctagonTriInds = rightOctagonTriInds;
malloid.rightOctagonVertInds = rightOctagonVertInds;




end















