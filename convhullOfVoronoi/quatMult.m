function z = quatMult(q1, q2)
s = q1(4);
v = q1(1:3);

t = q2(4);
w = q2(1:3);

mag = s*t - v*w';
dir = s*w + t*v + cross(v,w);

z = [dir mag];

end