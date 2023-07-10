


function h=visualize
    
    rands = rand(1000,3);
    close;
    hold on;
    
    r=[0,0,0];
    show(r, rands);
    
    
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
        if b1==source
            lr = source.Value;
        elseif b2 == source
            alt = source.Value;
        elseif b3 == source
            mag = source.Value;
        end
        
        
        alt
        lr
        r = mag*[cos(alt)*cos(lr) sin(alt)*cos(lr) sin(lr)]
        show(r, rands)
    end

end



function show(r, rands)
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

    % octahedron
    %{
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

    octb = [
        SO3([1,0,0], pi/4,rb);
        SO3([0,1,0], pi/4,rb);
        SO3([0,0,1], pi/4,rb);
        SO3([-1,0,0], pi/4,rb);
        SO3([0,-1,0], pi/4,rb);
        SO3([0,0,-1], pi/4,rb);
        
        SO3([-1,1,1], pi/3,rb);
        SO3([1,-1,1], pi/3,rb);
        SO3([1,1,-1], pi/3,rb);
        SO3([1,-1,-1], pi/3,rb);
        SO3([-1,-1,1], pi/3,rb);
        SO3([-1,1,-1], pi/3,rb);
        SO3([-1,-1,-1], pi/3,rb);
        SO3([1,1,1], pi/3,rb);
    ];
    %}
    %inds = [1 2 4 5 1 3 5 6 2 3 4 6 1 ];%8 5 8 3 11 4 11 5 10 1 10 6 9 2 9 1 2 12 4 12 6 13 4 13 5 3 14 1 14 2 7 3 7 4];
    %octmesh = oct(inds,:);
    %octmeshb = octb(inds,:);

    hold off; 
    scatter3(r(1),r(2),r(3), 50, [0,0,0]); 
    hold on; 
    pbaspect([1 1 1]); %axis off;
    plot3([0,r(1)],[0,r(2)],[0,r(3)],'color','k')
    %alpha(1);
    %plot3(octmesh(:,1),octmesh(:,2),octmesh(:,3),'color','r')
    %plot3(octmeshb(:,1),octmeshb(:,2),octmeshb(:,3),'color','k')
    
    %alpha(.1)
    
    
    for i = [1:size(axes,1)]
        plot3([axesb(i,1);axes(i,1)],[axesb(i,2);axes(i,2)],[axesb(i,3);axes(i,3)],'color','k');
    end
    
    for i = [1:size(corners,1)]
        plot3([cornersb(i,1);corners(i,1)],[cornersb(i,2);corners(i,2)],[cornersb(i,3);corners(i,3)],'color','k');
    end
    for i = [1:size(edges,1)]
        plot3([edgesb(i,1);edges(i,1)],[edgesb(i,2);edges(i,2)],[edgesb(i,3);edges(i,3)],'color','k');
    end
    %{
    for i = [1:size(oct,1)]
        plot3([oct(i,1);octb(i,1)],[oct(i,2);octb(i,2)],[oct(i,3);octb(i,3)],'color','k');
    end
    %}
    
    
    
    %scatter3(oct(:,1), oct(:,2), oct(:,3),30,[0,0,0]);
    a = scatter3(axes(:,1), axes(:,2), axes(:,3),30,[1,0,0]);
    c = scatter3(corners(:,1), corners(:,2), corners(:,3),30,[0,.5,0]);
    ident = scatter3(id(:,1), id(:,2), id(:,3),30,[0,0,0]);
    edg = scatter3(edges(:,1), edges(:,2), edges(:,3),30,[0,0,1]);
    
    
    %scatter3(axes(:,1), axes(:,2), axes(:,3),30,[1,0,0]);
    %scatter3(corners(:,1), corners(:,2), corners(:,3),30,[0,.7,0]);
    %scatter3(id(:,1), id(:,2), id(:,3),30,[0,0,0]);
    %scatter3(edges(:,1), edges(:,2), edges(:,3),30,[0,0,1]);
    
    
    plot3([-pi,pi],[0,0],[0,0],'color','r');
    plot3([0,0],[-pi,pi],[0,0],'color','b');
    plot3([0,0],[0,0],[-pi,pi],'color','g');
    
    [x,y,z]=sphere;
    Z = surf(pi*x,pi*y,pi*z,'edgecolor','none');
    alpha(Z,.05)
    
    rands = pi*rands./sqrt(sum(rands.^2,2));
    legend([a,c,ident,edg],{'axes','corner','identity','edge'});
    
    %h=scatter3(rands(:,1),rands(:,2),rands(:,3),1,[.8,.8,.8]);
    
end

          
          



%{
id = quaternion([1,1,1],0);
axes = [
quaternion([1,0,0], pi/2);
quaternion([1,0,0], pi);
quaternion([1,0,0], -pi/2);
quaternion([1,0,0], -pi);

quaternion([0,1,0], pi/2);
quaternion([0,1,0], pi);
quaternion([0,1,0], -pi/2);
quaternion([0,1,0], -pi);

quaternion([0,0,1], pi/2);
quaternion([0,0,1], pi);
quaternion([0,1,0], -pi/2);
quaternion([0,1,0], -pi);
];
edges = [
    quaternion([1,1,0],pi);
    quaternion([1,0,1],pi);
    quaternion([0,1,1],pi);
    
    quaternion([1,-1,0],pi);
    quaternion([1,0,-1],pi);
    quaternion([0,1,-1],pi);
    
    quaternion([-1,-1,0],pi);
    quaternion([-1,0,-1],pi);
    quaternion([0,-1,-1],pi);
    
    quaternion([-1,1,0],pi);
    quaternion([-1,0,1],pi);
    quaternion([0,-1,1],pi);
    ];
corners = [
quaternion([1,1,1],2*pi/3);
quaternion([1,1,1],-2*pi/3);
quaternion([1,1,-1],2*pi/3);
quaternion([1,1,-1],-2*pi/3);
quaternion([1,-1,1],2*pi/3);
quaternion([1,-1,1],-2*pi/3);
quaternion([-1,1,1],2*pi/3);
quaternion([-1,1,1],-2*pi/3);
quaternion([-1,-1,1],2*pi/3);
quaternion([-1,-1,1],-2*pi/3);
quaternion([1,-1,-1],2*pi/3);
quaternion([1,-1,-1],-2*pi/3);
quaternion([-1,1,-1],2*pi/3);
quaternion([-1,1,-1],-2*pi/3);
quaternion([-1,-1,-1],2*pi/3);
quaternion([-1,-1,-1],-2*pi/3);
];
%}