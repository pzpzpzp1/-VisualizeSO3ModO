% vertices in mesh1 and mesh2 that are in the exact same location are
% identified
function m3 = CombineMesh(mesh1,mesh2,Visualize,digits)
Visualize = 0;
if(nargin<4)
    digits = 5;
end

V1 = size(mesh1.V2P,1);
V2 = size(mesh2.V2P,1);
combinedV = [mesh1.V2P; mesh2.V2P];
[uniqueV, ia, ib] = unique(fix(combinedV*(10^digits)),'rows');
repeatInds = find(accumarray(ib,1)>1);

matchedInds = [];
for i = 1:numel(repeatInds)
    matchedInds = [matchedInds; find(ib(1:V1)==repeatInds(i)) find(ib(V1+1:end)==repeatInds(i))];
end

if Visualize
    figure; hold on; rotate3d on; scatter3(mesh1.V2P(:,1),mesh1.V2P(:,2),mesh1.V2P(:,3),'r');
    scatter3(mesh2.V2P(:,1),mesh2.V2P(:,2),mesh2.V2P(:,3),'r');
    scatter3(mesh1.V2P(matchedInds(:,1),1),mesh1.V2P(matchedInds(:,1),2),mesh1.V2P(matchedInds(:,1),3),10,'b');
    scatter3(mesh2.V2P(matchedInds(:,2),1),mesh2.V2P(matchedInds(:,2),2),mesh2.V2P(matchedInds(:,2),3),20,'k');
end

combinedV = [mesh1.V2P; mesh2.V2P()];
%combinedV(matchedInds(:,2)+V1,:)=[];

m3.V2P = combinedV;
m3.T2V = [mesh1.T2V; mesh2.T2V+V1];

matchedInds(:,2) = matchedInds(:,2) + V1;
Vs = 1:(V1+V2);
Vs(matchedInds(:,2)) = matchedInds(:,1);

m3.T2V = Vs(m3.T2V);

%% joint mesh made. time to compress. remove verts
verts2rem = matchedInds(:,2);

Vs = 1:(V1+V2);
Vs(matchedInds(:,2))=[];
V2s = 1:numel(Vs);
V3s = 1:numel(V2s);
V3s(Vs)=V2s;

m3.V2P(matchedInds(:,2),:)=[];
m3.T2V = V3s(m3.T2V);

if(Visualize)
    for i = 1:size(m3.T2V,1)
        tri = m3.T2V(i,:);
        ptc = m3.V2P(tri,:);
        ptc = patch(ptc(:,1),ptc(:,2),ptc(:,3),'green'); alpha(ptc,.3);
    end
end

end