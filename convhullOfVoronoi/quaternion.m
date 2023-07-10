function q = quaternion(dir, mag)
dir=dir/norm(dir);
q=[dir*sin(mag/2), cos(mag/2)];
end