%% test scale
clear;
n = 4;
dims = 2:1:min(10,n);
dist = rand(n); dist(1:n+1:end)=0; dist = dist + dist';
for i = dims
    C{i} = cmdscale(dist,i);
    M{i} = mdscale(dist,i);
end

[ii,jj]=find(ones(n)); e2v = sortrows([ii,jj]); e2v = e2v(find(e2v(:,1)<e2v(:,2)),:);
for i = dims
    V = C{i};
    EDists = dist;
    embeddedEdgeDifference = (V(e2v(:,1),:)-V(e2v(:,2),:));
    embeddedEdgeDistances = sqrt(sum(embeddedEdgeDifference.*embeddedEdgeDifference,2));
    idx = sub2ind(size(EDists), e2v(:,1), e2v(:,2));
    cscore = norm(EDists(idx)-embeddedEdgeDistances);
    
    V = M{i};
    embeddedEdgeDifference = (V(e2v(:,1),:)-V(e2v(:,2),:));
    embeddedEdgeDistances = sqrt(sum(embeddedEdgeDifference.*embeddedEdgeDifference,2));
    idx = sub2ind(size(EDists), e2v(:,1), e2v(:,2));
    mscore = norm(EDists(idx)-embeddedEdgeDistances);
    
    cscores(i)=cscore;
    mscores(i)=mscore;
end
cscores(dims)
mscores(dims)
    
