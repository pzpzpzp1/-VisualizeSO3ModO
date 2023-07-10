
function h=visualize
    
    rands = (rand(100000,3)-.5)*pi/2;
    
    rb = [0,0,0];
    idb = SO3([1,1,1],0,rb);
    axesb = [
        SO3([1,0,0], pi/2,rb);
        SO3([1,0,0], pi,rb);
        SO3([1,0,0], -pi/2,rb);
        SO3([1,0,0], -pi,rb);

        SO3([0,1,0], pi/2,rb);
        SO3([0,1,0], pi,rb);
        SO3([0,1,0], -pi/2,rb);
        SO3([0,1,0], -pi,rb);

        SO3([0,0,1], pi/2,rb);
        SO3([0,0,1], pi,rb);
        SO3([0,0,1], -pi/2,rb);
        SO3([0,0,1], -pi,rb);
    ];
    edgesb = [
        SO3([1,1,0],pi,rb);
        SO3([1,0,1],pi,rb);
        SO3([0,1,1],pi,rb);

        SO3([1,-1,0],pi,rb);
        SO3([1,0,-1],pi,rb);
        SO3([0,1,-1],pi,rb);

        SO3([-1,-1,0],pi,rb);
        SO3([-1,0,-1],pi,rb);
        SO3([0,-1,-1],pi,rb);

        SO3([-1,1,0],pi,rb);
        SO3([-1,0,1],pi,rb);
        SO3([0,-1,1],pi,rb);                         
    ];
    cornersb = [
        SO3([1,1,1],2*pi/3,rb);
        SO3([1,1,1],-2*pi/3,rb);
        SO3([1,1,-1],2*pi/3,rb);
        SO3([1,1,-1],-2*pi/3,rb);
        SO3([1,-1,1],2*pi/3,rb);
        SO3([1,-1,1],-2*pi/3,rb);
        SO3([-1,1,1],2*pi/3,rb);
        SO3([-1,1,1],-2*pi/3,rb);
        SO3([-1,-1,1],2*pi/3,rb);
        SO3([-1,-1,1],-2*pi/3,rb);
        SO3([1,-1,-1],2*pi/3,rb);
        SO3([1,-1,-1],-2*pi/3,rb);
        SO3([-1,1,-1],2*pi/3,rb);
        SO3([-1,1,-1],-2*pi/3,rb);
        SO3([-1,-1,-1],2*pi/3,rb);
        SO3([-1,-1,-1],-2*pi/3,rb);
    ];
    
    all = [idb;axesb;cornersb;edgesb];
    
    allu = unique([idb;axesb;cornersb;edgesb],'rows');
    %scatter3(allu(:,1),allu(:,2),allu(:,3),30,[1,0,0],'filled');
    [V,C] = voronoin(allu,{'Qbb'});
    x = V(:,1);y = V(:,2);z = V(:,3);
    
    
    for celli = 1:numel(C)
        celli;
        
        tcell = C{celli}; tcell(tcell==1)=tcell(2);
        %hold off; 
        scatter3(0,0,0,30,[0,0,0],'filled'); 
        hold on;
        scatter3(x(tcell),y(tcell),z(tcell),30,[1,0,0]);
        
        triangulationZ = convhull(x(tcell),y(tcell),z(tcell));
        surfTriang([x(tcell),y(tcell),z(tcell)], triangulationZ);
    end
    
    
    Stot = [];
    for iter=1:size(all,1)
        T = (rands-all(iter,:));
        S = sum(T.*T,2);
        Stot = [Stot,S];
    end
    rands = rands(min(Stot')'==Stot(:,1),:);
    triangulationZ = convhull(rands(:,1),rands(:,2),rands(:,3));
    rands = rands(unique(triangulationZ),:);
    triangulationZ = convhull(rands(:,1),rands(:,2),rands(:,3));
    
    
    
    close;
    hold on;
    
    r=[0,0,0];
    show(r, rands, triangulationZ);
    
    
    b1 = uicontrol('Style','slider','Position',[80,60,419,20], 'min',0, 'max',2*pi);
    b2 = uicontrol('Style','slider','Position',[80,80,419,20], 'min',0, 'max',pi);
    b3 = uicontrol('Style','slider','Position',[80,100,419,20], 'min',0, 'max',pi*2/3);
    b1.Callback=@temp;
    b2.Callback=@temp;
    b3.Callback=@temp;
    
    alt = 0;
    mag = 0;
    lr = 0;
    function temp(source,event)
        if b1==source
            lr = source.Value;
        elseif b2 == source
            alt = source.Value;
        elseif b3 == source
            mag = source.Value;
        end
        
        
        r = mag*[cos(alt)*cos(lr) sin(alt)*cos(lr) sin(lr)];
        r = mag*[1,1,1]/1.7321;
        show(r, rands, triangulationZ)
    end

end

function showOcts(r,color)
    oct = [
        SO3([1,0,0], pi/4,r);
        SO3([0,1,0], pi/4,r);
        SO3([0,0,1], pi/4,r);
        SO3([-1,0,0], pi/4,r);
        SO3([0,-1,0], pi/4,r);
        SO3([0,0,-1], pi/4,r);
        
        
        SO3([-1,1,1], pi/3,r);
        SO3([1,-1,1], pi/3,r);
        SO3([1,1,-1], pi/3,r);
        SO3([1,-1,-1], pi/3,r);
        SO3([-1,-1,1], pi/3,r);
        SO3([-1,1,-1], pi/3,r);
        SO3([-1,-1,-1], pi/3,r);
        SO3([1,1,1], pi/3,r);
        
    ];

    inds = [1 2 4 5 1 3 5 6 2 3 4 6 1     8 5 8 3 11 4 11 5 10 1 10 6 9 2 9 1 2 12 4 12 6 13 4 13 5 3 14 1 14 2 7 3 7 4];
    %inds = [1 2 4 5 1 3 5 6 2 3 4 6 1];%   8 5 8 3 11 4 11 5 10 1 10 6 9 2 9 1 2 12 4 12 6 13 4 13 5 3 14 1 14 2 7 3 7 4];
    
    octmesh = oct(inds,:);
    C = [0.5000 1.0000 1.0000 0.5000;
     1.0000 0.5000 0.5000 0.1667;
     0.3330 0.3330 0.5000 0.5000];
    plot3(octmesh(:,1),octmesh(:,2),octmesh(:,3),'color',color);
    
    for iter = 1:size(oct,1)
        scatter3(oct(iter,1),oct(iter,2),oct(iter,3),30,[1,0,0],'filled');
    end
    
end

function show(r, rands, triangulationZ)
    id = SO3([1,1,1],0,r);
    axes = [
        SO3([1,0,0], pi/2,r);
        SO3([1,0,0], pi,r);
        SO3([1,0,0], -pi/2,r);
        SO3([1,0,0], -pi,r);

        SO3([0,1,0], pi/2,r);
        SO3([0,1,0], pi,r);
        SO3([0,1,0], -pi/2,r);
        SO3([0,1,0], -pi,r);

        SO3([0,0,1], pi/2,r);
        SO3([0,0,1], pi,r);
        SO3([0,0,1], -pi/2,r);
        SO3([0,0,1], -pi,r);
    ];
    edges = [
        SO3([1,1,0],pi,r);
        SO3([1,0,1],pi,r);
        SO3([0,1,1],pi,r);

        SO3([1,-1,0],pi,r);
        SO3([1,0,-1],pi,r);
        SO3([0,1,-1],pi,r);

        SO3([-1,-1,0],pi,r);
        SO3([-1,0,-1],pi,r);
        SO3([0,-1,-1],pi,r);

        SO3([-1,1,0],pi,r);
        SO3([-1,0,1],pi,r);
        SO3([0,-1,1],pi,r);                         
    ];
    corners = [
        SO3([1,1,1],2*pi/3,r);
        SO3([1,1,1],-2*pi/3,r);
        SO3([1,1,-1],2*pi/3,r);
        SO3([1,1,-1],-2*pi/3,r);
        SO3([1,-1,1],2*pi/3,r);
        SO3([1,-1,1],-2*pi/3,r);
        SO3([-1,1,1],2*pi/3,r);
        SO3([-1,1,1],-2*pi/3,r);
        SO3([-1,-1,1],2*pi/3,r);
        SO3([-1,-1,1],-2*pi/3,r);
        SO3([1,-1,-1],2*pi/3,r);
        SO3([1,-1,-1],-2*pi/3,r);
        SO3([-1,1,-1],2*pi/3,r);
        SO3([-1,1,-1],-2*pi/3,r);
        SO3([-1,-1,-1],2*pi/3,r);
        SO3([-1,-1,-1],-2*pi/3,r);
    ];

    


    rb = [0,0,0];
    idb = SO3([1,1,1],0,rb);
    axesb = [
        SO3([1,0,0], pi/2,rb);
        SO3([1,0,0], pi,rb);
        SO3([1,0,0], -pi/2,rb);
        SO3([1,0,0], -pi,rb);

        SO3([0,1,0], pi/2,rb);
        SO3([0,1,0], pi,rb);
        SO3([0,1,0], -pi/2,rb);
        SO3([0,1,0], -pi,rb);

        SO3([0,0,1], pi/2,rb);
        SO3([0,0,1], pi,rb);
        SO3([0,0,1], -pi/2,rb);
        SO3([0,0,1], -pi,rb);
    ];
    edgesb = [
        SO3([1,1,0],pi,rb);
        SO3([1,0,1],pi,rb);
        SO3([0,1,1],pi,rb);

        SO3([1,-1,0],pi,rb);
        SO3([1,0,-1],pi,rb);
        SO3([0,1,-1],pi,rb);

        SO3([-1,-1,0],pi,rb);
        SO3([-1,0,-1],pi,rb);
        SO3([0,-1,-1],pi,rb);

        SO3([-1,1,0],pi,rb);
        SO3([-1,0,1],pi,rb);
        SO3([0,-1,1],pi,rb);                         
    ];
    cornersb = [
        SO3([1,1,1],2*pi/3,rb);
        SO3([1,1,1],-2*pi/3,rb);
        SO3([1,1,-1],2*pi/3,rb);
        SO3([1,1,-1],-2*pi/3,rb);
        SO3([1,-1,1],2*pi/3,rb);
        SO3([1,-1,1],-2*pi/3,rb);
        SO3([-1,1,1],2*pi/3,rb);
        SO3([-1,1,1],-2*pi/3,rb);
        SO3([-1,-1,1],2*pi/3,rb);
        SO3([-1,-1,1],-2*pi/3,rb);
        SO3([1,-1,-1],2*pi/3,rb);
        SO3([1,-1,-1],-2*pi/3,rb);
        SO3([-1,1,-1],2*pi/3,rb);
        SO3([-1,1,-1],-2*pi/3,rb);
        SO3([-1,-1,-1],2*pi/3,rb);
        SO3([-1,-1,-1],-2*pi/3,rb);
    ];


    
        hold off; 
        %subplot(1,3,jj);
        
        scatter3(r(1),r(2),r(3), 50, [0,0,0]); 
        hold on;
        pbaspect([1 1 1]); %axis off;
        plot3([0,r(1)],[0,r(2)],[0,r(3)],'color','k')

        showOcts(idb(1,:),'k');

        scatter3(axesb(:,1), axesb(:,2), axesb(:,3),10,[1,0,0]);
        scatter3(cornersb(:,1), cornersb(:,2), cornersb(:,3),10,[1,0,0]);
        scatter3(idb(:,1), idb(:,2), idb(:,3),10,[1,0,0]);
        scatter3(edgesb(:,1), edgesb(:,2), edgesb(:,3),10,[1,0,0]);

        plot3([-pi,pi],[0,0],[0,0],'color','r');
        plot3([0,0],[-pi,pi],[0,0],'color','b');
        plot3([0,0],[0,0],[-pi,pi],'color','g');

        [x,y,z]=sphere;
        Z = surf(pi*x,pi*y,pi*z,'edgecolor','none');
        alpha(Z,.05)
        
        beg=16; ed=21; 
        Z2 = surf(pi*abs(x(beg:ed,beg:ed)),pi*y(beg:ed,beg:ed),pi*z(beg:ed,beg:ed),'edgecolor','none');
        alpha(Z2, .6)
        
        %grid off; axis off;
        %rands = pi*rands./sqrt(sum(rands.^2,2));
        %h=scatter3(rands(:,1),rands(:,2),rands(:,3),1,[.8,.8,.8]);
        surfTriang(rands, triangulationZ);

end

function surfTriang(points, triang)
    hold on;
    x = points(:,1); y = points(:,2); z= points(:,3);
    ptc = patch(x(triang'),y(triang'),z(triang'),'green');
    alpha(ptc, .05);
    %hold on;
    %for p = 1:size(triang,1)
    %    surfT(points(triang(p,1),:),points(triang(p,2),:),points(triang(p,3),:));
    %end
end

function surfT(p1,p2,p3)
    p4=(p1+p2)/2;

    X = [p1(1),p2(1);p3(1),p4(1)];
    Y = [p1(2),p2(2);p3(2),p4(2)];
    Z = [p1(3),p2(3);p3(3),p4(3)];
    
    %patch(p1(1),p2(1),p3(1))
    
    surf(X,Y,Z);
    
end
          
