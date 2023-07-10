% empirical sampling of bisecting rotation between identity and qa
function z = Bisector(aa,n)
    if(nargin==0)
        n=100000;
        aa = [0,0,pi/2];
    end
    
    
    qb = axang2quat([aa/norm(aa) norm(aa)/2]);
    dirs = randn(n,3)*100; dirs = dirs ./ sqrt(sum(dirs.*dirs,2));
    mags = mod(randn(n,1)*100000,2*pi);
    qs = axang2quat([mags dirs]);
    qconjs = quatmultiply(quatmultiply(quatinv(qs),repmat(qb,n,1)),qs);
    conjugatedAAs = quat2axang(qconjs);
    conjugatedAAs = conjugatedAAs(:,1:3).*conjugatedAAs(:,4);
    z = conjugatedAAs;
    
end