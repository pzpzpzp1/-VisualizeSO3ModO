% without mod!!! more interesting than default quat2axang
function z = quatToAxisAngle(q)

qijk = q(:,2:4);
qr = q(:,1);

norms = sqrt(sum(qijk.*qijk,2));
dir = qijk./norms;
mag = 2*atan2(norms, qr);

z = [dir mag];
z(not(z==z))=0;
end