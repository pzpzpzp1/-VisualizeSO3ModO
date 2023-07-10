
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
        pause
    end
    
    
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
          
