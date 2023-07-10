function z = quatToAxisAngle(q)

qijk = q(1:3);
qr = q(4);

dir = qijk/norm(qijk);
mag = 2*atan2(norm(qijk), qr);

z = [dir mag];
z(not(z==z))=0;
end