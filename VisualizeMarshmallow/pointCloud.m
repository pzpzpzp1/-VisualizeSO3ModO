addpath('Marshmallow');
addpath('..\..\SHFramesInVolume');
addpath('../../../jsolomon/octahedral_frames/code/harmonics');
addpath('../../../jsolomon/octahedral_frames/code/frames');
addpath('..\..\..\jsolomon\octahedral_frames\code\external\easyspin-5.0.11\easyspin');
surfaceResolution = 6;
idxs = [1 2 3];
mdim = 9;

% axis aligned parameters
v = [1 1 1];
n = 40;

malloid = Marshmallow(surfaceResolution,0);
nSurfacePoints = size(malloid.V2P,1)
innerResolution = ceil((sqrt(size(malloid.V2P,1)/6))^3)
points = [];
while (size(points,1) < innerResolution)
    pointsToAdd = (rand(innerResolution,3)-.5)*pi/2;
    [z,d,inds] = equivalenceInSO3ModO(pointsToAdd, zeros(innerResolution,3));
    points = [points; pointsToAdd(find(inds==1),:)];
end; points = points(1:innerResolution,:);
points = [malloid.V2P; points];

% points showing rotations aligned with v
AA = quatToAxisAngle(axisRestrictedQuats(v, n));
aa = AA(:,1:3).*AA(:,4);
bb = aa;
[zs,ds,inds,aa] = equivalenceInSO3ModO(aa, repmat([0 0 0],size(aa,1),1));
[zs,ds,inds,aa] = equivalenceInSO3ModO(aa, repmat([0 0 0],size(aa,1),1));
assert(all(inds));
size(points,1)
points = [points; aa];
AAinds = (size(points,1)-size(aa,1))+1:size(points,1);


% compute pairwise dists
nV1 = size(points,1);
i1 = repmat(1:nV1,1,nV1)';
i2 = repmat(1:nV1,nV1,1); i2 = i2(:);
removeInds = find(i1<=i2);
i1(removeInds)=[]; i2(removeInds)=[];
[~,ds,~] = equivalenceInSO3ModO(points(i1,:),points(i2,:));
dist = zeros(nV1,nV1);
dist(sub2ind([nV1,nV1], i1, i2)) = ds;
dist = dist + dist';

%% SPH comparison
% compute sph score
Vsph = zeros(size(points,1),9);
for i = 1:size(points,1)
    Vsph(i,:) = .5*rotm2sph(axang2rotm([points(i,:) norm(points(i,:))]));
end
dV = Vsph(i1,:)-Vsph(i2,:);
dVnorm = sqrt(sum(dV.*dV,2));
% score1 is 301.5, score2 is .2044. at optimal scaling.
% scores after embedding. 188 .1374
% originalD = dist(sub2ind([nV1,nV1], i1, i2));
% sphscore1 = [];
% sphscore2 = [];
% for scale = 1:200
%     scl = scale/100;
%     sphscore1(scale) = norm(dVnorm*scl-originalD);
%     sphscore2(scale) = mean(abs(dVnorm*scl-originalD))/mean([dVnorm*scl;originalD]);
% end
% figure; hold on; plot(sphscore1,'r'); title('sphscore1');
% figure; hold on; plot(sphscore2,'g'); title('sphscore2');

'computing embedding';
tic
    % V = cmdscale(dist);

    % dist3 = dist;
    % dist3(dist3>.2)=nan;
    % V0 = cmdscale(dist,mdim);
    % V = mdscale(dist3,mdim,'start',V0);

    dist2 = dist;
    ii = find(dist2==0); dist2(ii) = .0001*rand(size(dist2(ii)));
    dist2(1:(size(dist2,1)+1):end)=0;
    dist2 = (dist2 + dist2') /2;
    V = mdscale(dist2,mdim,'start',Vsph+.0001*rand(size(Vsph)));
    %V = mdscale(dist2,mdim);
    
%     V = Vsph
toc



score1 = []; score2 = [];
for idim = 1:mdim;
    dV = V(i1,1:idim )-V(i2,1:idim );
    dVnorm = sqrt(sum(dV.*dV,2));
    originalD = dist(sub2ind([nV1,nV1], i1, i2));
    score1(idim ) = norm(dVnorm-originalD);
    score2(idim ) = mean(abs(dVnorm-originalD))/mean([dVnorm;originalD]);
end
figure; hold on; plot(score1,'r'); title('score1');
figure; hold on; plot(score2,'g'); title('score2');

% compute colormaps
Cra = sqrt(sum(points.*points,2));
Caz = atan2(points(:,1),points(:,2));
rescaledCra = (Cra - min(Cra))./max(Cra - min(Cra));
rescaledCaz = (Caz - min(Caz))./max(Caz - min(Caz));

% show embedding
f = figure; hold all; axis equal; rotate3d on; title('colored radially initial');
scatter3(points(:,1),points(:,2),points(:,3),5,rescaledCra);
PlotMallow(malloid, 0, f, 'green',.05,'none');

f = figure; hold all; axis equal; rotate3d on; title('colored azimuthally initial');
scatter3(points(:,1),points(:,2),points(:,3),5,rescaledCaz);
PlotMallow(malloid, 0, f, 'green',.05,'none');

figure; hold all; axis equal; rotate3d on; title('colored radially embedded cmds');
scatter3(V(:,idxs(1)),V(:,idxs(2)),V(:,idxs(3)),5,rescaledCra);

figure; hold all; axis equal; rotate3d on; title('colored azimuthally embedded cmds');
scatter3(V(:,idxs(1)),V(:,idxs(2)),V(:,idxs(3)),5,rescaledCaz);

figure; hold all; axis equal; rotate3d on; title('show edges corners and faces colorcoded embedded cmds');
scatter3(V(:,idxs(1)),V(:,idxs(2)),V(:,idxs(3)),5,'k.');
scatter3(V(malloid.faceInds,idxs(1)),V(malloid.faceInds,idxs(2)),V(malloid.faceInds,idxs(3)),'g');
scatter3(V(malloid.edgeInds,idxs(1)),V(malloid.edgeInds,idxs(2)),V(malloid.edgeInds,idxs(3)),'b');
scatter3(V(malloid.cornerInds,idxs(1)),V(malloid.cornerInds,idxs(2)),V(malloid.cornerInds,idxs(3)),500,'k.');

%% animated version
C = sin(rescaledCaz*2*pi);
V0 = points;
V1 = V;
figure; hold all; axis equal; rotate3d on; title('colored RA embedded cmds animated');
c = 1;
dt = .01;
while(c < 100); c = c+1;
    for t = [[0:dt:1] fliplr([0:dt:1])]
        cla;
        Vt = V0*(1-t)+V1(:,idxs)*t;
        
        %% draw verts
        scatter3(Vt(:,1),Vt(:,2),Vt(:,3),5,C);
        %scatter3(Vt(malloid.edgeInds,1),Vt(malloid.edgeInds,2),Vt(malloid.edgeInds,3),'b');
        scatter3(Vt(malloid.cornerInds,1),Vt(malloid.cornerInds,2),Vt(malloid.cornerInds,3),500,'k.');
        %scatter3(V(malloid.faceInds,idxs(1)),V(malloid.faceInds,idxs(2)),V(malloid.faceInds,idxs(3)),'g');
        
        %% draw lines and faces
        patch('Faces',malloid.SE2V,'Vertices',Vt);
        %ptc = patch('Faces',malloid.T2V,'Vertices',Vt,'FaceColor','green'); alpha(ptc,.1);
        ptc = patch('Faces',malloid.T2V(malloid.rightOctagonTriInds,:),'Vertices',Vt,'FaceColor','green','EdgeColor','k'); alpha(ptc,.1);

        %% draw axis aligned rotations
        scatter3(Vt(AAinds,1),Vt(AAinds,2),Vt(AAinds,3),'r');
        
        drawnow;
        if(t==1 || t==0)
            pause
        end
    end
end