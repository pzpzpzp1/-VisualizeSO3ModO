function [Vout,Qout] = shiftPoints(V,q)
    assert(size(V,2)==3);
    assert(numel(q)==4);
    
    n = size(V,1);
    norms = sqrt(sum(V.*V,2));
    V(find(norms==0),1) = V(find(norms==0),1) + 1;
    
    qs = axang2quat([V norms]);
    
    q2s = quatmultiply(qs, repmat(q(:)',n,1));
    
    %% these two are different. quat2axang is modded by pi in magnitude.
    %AAout = quat2axang(q2s);
    AAout = quatToAxisAngle(q2s);
    
    Vout = AAout(:,1:3).*AAout(:,4);
    Qout = q2s;
end