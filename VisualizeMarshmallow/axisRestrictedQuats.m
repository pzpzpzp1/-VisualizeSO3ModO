function quats = axisRestrictedQuats(v, n)
    if nargin==0
        v = [1 2 1]';
        n = 10;
    end

    v = v(:);
    assert(numel(v)==3);
    assert(numel(n)==1);
    assert(n > 1);
    v = v/norm(v);
    ang = linspace(0,pi/2,n+1)';
    quats = [];
    
    start = [0 0 1]';
    Rq = rotm2quat(FromStartToEnd(start, v));
    Rz = axang2quat([repmat(start',n,1) ang(1:end-1)]);
    quats = [quats; quatmultiply(repmat(Rq,n,1),Rz)];
    
    start = [0 1 0]';
    Rq = rotm2quat(FromStartToEnd(start, v));
    Ry = axang2quat([repmat(start',n,1) ang(1:end-1)]);
    quats = [quats; quatmultiply(repmat(Rq,n,1),Ry)];

    start = [1 0 0]';
    Rq = rotm2quat(FromStartToEnd(start, v));
    Rx = axang2quat([repmat(start',n,1) ang(1:end-1)]);
    quats = [quats; quatmultiply(repmat(Rq,n,1),Rx)];
    
    start = [0 0 -1]';
    Rq = rotm2quat(FromStartToEnd(start, v));
    Rmz = axang2quat([repmat(start',n,1) ang(1:end-1)]);
    quats = [quats; quatmultiply(repmat(Rq,n,1),Rmz)];
    
    start = [0 -1 0]';
    Rq = rotm2quat(FromStartToEnd(start, v));
    Rmy = axang2quat([repmat(start',n,1) ang(1:end-1)]);
    quats = [quats; quatmultiply(repmat(Rq,n,1),Rmy)];
    
    start = [-1 0 0]';
    Rq = rotm2quat(FromStartToEnd(start, v));
    Rmx = axang2quat([repmat(start',n,1) ang(1:end-1)]);
    quats = [quats; quatmultiply(repmat(Rq,n,1),Rmx)];
    
end

function R = FromStartToEnd(startVector, v)
    if(startVector'*v == -1)
        vs = null([0 1 0]); 
        vs = vs(:,1)';
        R = axang2rotm([vs pi]);
    else
        vec = cross(startVector', v');
        c = startVector' * v;
        vm = [0 -vec(3) vec(2); vec(3) 0 -vec(1); -vec(2) vec(1) 0];
        R = eye(3) + vm + vm*vm/(1+c);
    end
    assert(norm(R*startVector-v)<.0000001);
end


