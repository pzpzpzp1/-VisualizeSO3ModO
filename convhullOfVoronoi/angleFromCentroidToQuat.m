% takes two vectors in SO3 representing quaternions, returns angle
% (distance) between them.
function z = angleFromCentroidToQuat(centroidSO3, destSO3)
    
quatCentroid = quaternion(centroidSO3, norm(centroidSO3));
quatDest = quaternion(destSO3, norm(destSO3));

qx = quatdivide([quatDest(4) quatDest(1:3)],[quatCentroid(4) quatCentroid(1:3)]);

axang = quat2axang(qx); % this is in SO3!

z = axang(4);
end