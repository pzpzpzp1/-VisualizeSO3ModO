
function z = SO3(dir, mag, r)
dir=dir/norm(dir);
%q=[dir*mag];

rmag = norm(r);
r = r/rmag;
if(rmag==0)
    r=[1,0,0];
end
    
q = quatMult(quaternion(dir, mag),quaternion(r, rmag));
AA = quatToAxisAngle(q);
z = AA(1:3)*AA(4);
end