function SaveMalloidAsSmesh(fid, malloid)
    V2P = malloid.V2P;
    F2V = malloid.T2V;
    
    n = size(V2P,1);
    F = size(F2V,1);
    Ps = [(1:n)',V2P]';
    fprintf(fid,'%d %d 0 0\n',n,3);
    fprintf(fid,'%d %d %d %d\n',Ps(:));
    
    Fs = F2V';
    fprintf(fid,'%d 0\n',F);
    fprintf(fid,'3 %d %d %d\n',Fs(:));
    
end