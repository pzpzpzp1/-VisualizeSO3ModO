addpath('Marshmallow');
addpath('..\..\SHFramesInVolume');
addpath('../../../jsolomon/octahedral_frames/code/harmonics');
addpath('../../../jsolomon/octahedral_frames/code/frames');
addpath('..\..\..\jsolomon\octahedral_frames\code\external\easyspin-5.0.11\easyspin');

n=100000;
record = zeros(1,n);
for i = 1:n
    % rotation
    AR = eye(3);
    [BR,~] = qr(rand(3));
    [CR,~] = qr(rand(3));
    
    AAA = rotm2axang(AR);
    AAB = rotm2axang(BR);
    AAC = rotm2axang(CR);
    
    aaA = AAA(1:3)*AAA(4);
    aaB = AAB(1:3)*AAB(4);
    aaC = AAC(1:3)*AAC(4);
    
    [~,~,~,aaA] = equivalenceInSO3ModO(aaA, [0,0,0]);
    [~,~,indA,aaA] = equivalenceInSO3ModO(aaA, [0,0,0]);
    [~,~,~,aaB] = equivalenceInSO3ModO(aaB, [0,0,0]);
    [~,~,indB,aaB] = equivalenceInSO3ModO(aaB, [0,0,0]);
    [~,~,~,aaC] = equivalenceInSO3ModO(aaC, [0,0,0]);
    [~,~,indC,aaC] = equivalenceInSO3ModO(aaC, [0,0,0]);
    assert(all([indA indB indC]));
    
    [~,dAB] = equivalenceInSO3ModO(aaA, aaB);
    [~,dAC] = equivalenceInSO3ModO(aaA, aaC);
    
    sphA = rotm2sph(AR);
    sphB = rotm2sph(BR);
    sphC = rotm2sph(CR);
    
    dsphAB = norm(sphA-sphB);
    dsphAC = norm(sphA-sphC);
    
    record(i) = (dAB <= dAC) == (dsphAB <= dsphAC);
end

maintainsDistance = sum(record)/n

