function ComputeAbstractMalloids()
    assert(exist('cache')~=0)
    addpath('Marshmallow')
    addpath('AbstractMalloid')
    
    for dimension = 3:8
        for resolution = [20:10:40]
            dimension
            resolution
            
            fnameM = ['cache/malloidD' num2str(dimension) 'R' num2str(resolution) '.mat'];
            fnameMT = ['cache/malloidTetmeshesD' num2str(dimension) 'R' num2str(resolution) '.mat'];
            fnameMAT = ['cache/malloidAbstractTetmeshesD' num2str(dimension) 'R' num2str(resolution) '.mat'];

            if(exist(fnameM)==0)
                [malloid, malloidTetmesh, malloidAbstractTetmesh]= Marshmallow(resolution, 0, dimension);
            
                
                'saving'
                tic
                save(fnameM,'malloid','-v7.3');
                save(fnameMT,'malloidTetmesh','-v7.3');
                save(fnameMAT,'malloidAbstractTetmesh','-v7.3');
                toc
            end
        end
    end
end