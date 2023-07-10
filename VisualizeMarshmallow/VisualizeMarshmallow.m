
function h = VisualizeMarshmallows
    
    addpath('Marshmallow');
    res = 11;
    malloid = Marshmallow(res,0);

    mallowname = ['malloid' num2str(res) '.mat'];
    if(exist(mallowname))
        load(mallowname);
    else
        save(mallowname, 'malloid');
    end
    [id,axes,edges,corners,all] = hardcodedAA();
    
    f = figure; hold on; rotate3d on; axis equal; axis off; grid off;
    zbox = 2*pi; xlim([-zbox,zbox]); ylim([-zbox,zbox]); zlim([-zbox,zbox]);
    showAxis();
    showInds = ones(size(all,1),1)*0;
    % number of values per type. 1 24 24 48 id,axes,edges,corners
    showInds = [1, ...
      1,1,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,...
      1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,...
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
    
%     % Ball within pi
%     showInds = [1, ...
%       1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,...
%       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,...
%       1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
    
    % Ball on Pi
    showInds = [0, ...
      0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,...
      1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,...
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
    showCells(malloid, f, showInds, [1 0 0 0]);
    
    b1 = uicontrol('Style','slider','Position',[80,60,419,20], 'min',0, 'max',2*pi);
    b2 = uicontrol('Style','slider','Position',[80,80,419,20], 'min',0, 'max',pi);
    b3 = uicontrol('Style','slider','Position',[80,100,419,20], 'min',0, 'max',pi);
    b1.Callback=@temp;
    b2.Callback=@temp;
    b3.Callback=@temp;
    
    alt = 0;
    mag = 0;
    lr = 0;
    function temp(source,event)
        cla; showAxis(); xlim([-zbox,zbox]); ylim([-zbox,zbox]); zlim([-zbox,zbox]);
        
        if b1==source
            lr = source.Value;
        elseif b2 == source
            alt = source.Value;
        elseif b3 == source
            mag = source.Value;
        end
        
        r = mag*[cos(alt)*cos(lr) sin(alt)*cos(lr) sin(lr)];
        plot3([0,r(1)],[0,r(2)],[0,r(3)],'color','k','LineWidth',3)
        
        Q = axang2quat([r norm(r)]);
        if(norm(r)==0)
            Q = [1,0,0,0];
        end
        
        %% toggle this to make cells move w/w/o ui controls
        %showCells(malloid, f, showInds, [1 0 0 0]);
        showCells(malloid, f, showInds, Q);
        
        PlotMallow(malloid,Q,f,'red',-1);
    end
    
    
end

function showAxis()
    plot3([-pi,pi],[0,0],[0,0],'color','r');
    plot3([0,0],[-pi,pi],[0,0],'color','b');
    plot3([0,0],[0,0],[-pi,pi],'color','g');
    [x,y,z]=sphere;
    Z = surf(pi*x,pi*y,pi*z,'edgecolor','none'); alpha(Z,.05)
    Z = surf(2*pi*x,2*pi*y,2*pi*z,'edgecolor','none'); alpha(Z,.05)
end

function showCells(malloid, f, showInds, quat);
    [id,axes,corners,edges,all] = hardcodedAA();
    aacenters = all(find(showInds),:);
    if(numel(aacenters)==0); return; end;
    [V, Qout]  = shiftPoints(aacenters(:,1:3).*aacenters(:,4),quat);
    for i = 1:size(Qout,1)
        PlotMallow(malloid,Qout(i,:),f,0,-1);
    end
end

         

function [id,axes,edges,corners,all] = hardcodedAA()
    id = [1,1,1,0];
    axes = [
        [1,0,0,pi/2];
        [1,0,0,pi];
        [1,0,0,3/2*pi];
        [1,0,0,2*pi];
        
        [1,0,0,-pi/2];
        [1,0,0,-pi];
        [1,0,0,-3/2*pi];
        [1,0,0,-2*pi];

        [0,1,0,pi/2];
        [0,1,0,pi];
        [0,1,0,3/2*pi];
        [0,1,0,2*pi];
        
        [0,1,0,-pi/2];
        [0,1,0,-pi];
        [0,1,0,-3/2*pi];
        [0,1,0,-2*pi];

        [0,0,1,pi/2];
        [0,0,1,pi];
        [0,0,1,3/2*pi];
        [0,0,1,2*pi];
        
        [0,0,1,-pi/2];
        [0,0,1,-pi];
        [0,0,1,-3/2*pi];
        [0,0,1,-2*pi];
    ];
    edges = [
        [1,1,0,pi];
        [1,0,1,pi];
        [0,1,1,pi];

        [1,-1,0,pi];
        [1,0,-1,pi];
        [0,1,-1,pi];

        [-1,-1,0,pi];
        [-1,0,-1,pi];
        [0,-1,-1,pi];

        [-1,1,0,pi];
        [-1,0,1,pi];
        [0,-1,1,pi];
        
        % 2pi
        [1,1,0,2*pi];
        [1,0,1,2*pi];
        [0,1,1,2*pi];

        [1,-1,0,2*pi];
        [1,0,-1,2*pi];
        [0,1,-1,2*pi];

        [-1,-1,0,2*pi];
        [-1,0,-1,2*pi];
        [0,-1,-1,2*pi];

        [-1,1,0,2*pi];
        [-1,0,1,2*pi];
        [0,-1,1,2*pi];
    ];
    corners = [
        [1,1,1,2*pi/3];
        [1,1,1,-2*pi/3];
        [1,1,-1,2*pi/3];
        [1,1,-1,-2*pi/3];
        [1,-1,1,2*pi/3];
        [1,-1,1,-2*pi/3];
        [-1,1,1,2*pi/3];
        [-1,1,1,-2*pi/3];
        [-1,-1,1,2*pi/3];
        [-1,-1,1,-2*pi/3];
        [1,-1,-1,2*pi/3];
        [1,-1,-1,-2*pi/3];
        [-1,1,-1,2*pi/3];
        [-1,1,-1,-2*pi/3];
        [-1,-1,-1,2*pi/3];
        [-1,-1,-1,-2*pi/3];
        
        [1,1,1,4*pi/3];
        [1,1,1,-4*pi/3];
        [1,1,-1,4*pi/3];
        [1,1,-1,-4*pi/3];
        [1,-1,1,4*pi/3];
        [1,-1,1,-4*pi/3];
        [-1,1,1,4*pi/3];
        [-1,1,1,-4*pi/3];
        [-1,-1,1,4*pi/3];
        [-1,-1,1,-4*pi/3];
        [1,-1,-1,4*pi/3];
        [1,-1,-1,-4*pi/3];
        [-1,1,-1,4*pi/3];
        [-1,1,-1,-4*pi/3];
        [-1,-1,-1,4*pi/3];
        [-1,-1,-1,-4*pi/3];
        
        [1,1,1,6*pi/3];
        [1,1,1,-6*pi/3];
        [1,1,-1,6*pi/3];
        [1,1,-1,-6*pi/3];
        [1,-1,1,6*pi/3];
        [1,-1,1,-6*pi/3];
        [-1,1,1,6*pi/3];
        [-1,1,1,-6*pi/3];
        [-1,-1,1,6*pi/3];
        [-1,-1,1,-6*pi/3];
        [1,-1,-1,6*pi/3];
        [1,-1,-1,-6*pi/3];
        [-1,1,-1,6*pi/3];
        [-1,1,-1,-6*pi/3];
        [-1,-1,-1,6*pi/3];
        [-1,-1,-1,-6*pi/3];
    ];
    all = [id;axes;edges;corners];
    
    id(:,1:3) = id(:,1:3)./sqrt(sum(id(:,1:3).*id(:,1:3),2));
    axes(:,1:3) = axes(:,1:3)./sqrt(sum(axes(:,1:3).*axes(:,1:3),2));
    edges(:,1:3) = edges(:,1:3)./sqrt(sum(edges(:,1:3).*edges(:,1:3),2));
    corners(:,1:3) = corners(:,1:3)./sqrt(sum(corners(:,1:3).*corners(:,1:3),2));
    all(:,1:3) = all(:,1:3)./sqrt(sum(all(:,1:3).*all(:,1:3),2));
end
