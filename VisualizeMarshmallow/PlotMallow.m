function fig = PlotMallow(malloid, q, f, color, aparam, linestyle)
    if(f==0)
        fig = figure; hold on; rotate3d on; axis equal;
    else
        fig = figure(f);
    end

    if(q==0)
        q = [1 0 0 0];
    end
    
    if(color==0)
        color = 'green'
    end
    
    if(aparam == -1 || nargin < 5)
        aparam = .3;
    end
    
    if(nargin < 6)
        linestyle = '-';
    end

    malloid.V2P = shiftPoints(malloid.V2P, q);

    ptc = patch('Faces',malloid.T2V,'Vertices',malloid.V2P,'FaceColor',color,'LineStyle',linestyle); alpha(ptc,aparam);
end