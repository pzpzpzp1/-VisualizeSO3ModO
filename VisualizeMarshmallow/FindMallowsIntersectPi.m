%% some theory. not that much theory... a wedge on the face of the mallow is parameterized in quaternion space by 3 points in R4. Any partition of union and positivity interpolation of them yields the wedge face, and can be converted into AA space (R3) to make curved surfaces.
%% in trying to find where the surfaces intersect Pi, we just need the first value in R4 representation to be 0. Quaternion first element is cos(theta/2), and cos(pi/2)=0.
%% Similarly, in trying to find where the surfaces intersect 2Pi, we just need the last three values in R4 representation to be 0. Quaternion last three elements are scaled by sin(theta/2), and sin(pi)=0.

% BTW THE TRIANGLE FACES ARE SOLID!!! WHICH IS NEAT.

function [aaverts, solutionQEdges] = FindMallowsIntersectPi
    
    addpath('Marshmallow');
    lineres = 10;
    malloid2 = Marshmallow(2,0);
    malloid11 = Marshmallow(11,0);

    % number of values per type. 1 24 24 48 id,axes,edges,corners
    % ball on pi
    showInds = [0, ...
      0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,...
      1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,...
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
    showInds = [0, ...
      0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,...
      1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,...
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
%     showInds = [0, ...
%       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,...
%       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,...
%       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
    
    aaverts = [];
    [id,axes,corners,edges,all] = hardcodedAA();
    aacenters = all(find(showInds),:);
    solutionQEdges = [];
    
    f = figure; hold on; axis equal; rotate3d on;
    [x,y,z]=sphere;
    Z = surf(pi*x,pi*y,pi*z,'edgecolor','none'); alpha(Z,.8)
    %Z = surf(2*pi*x,2*pi*y,2*pi*z,'edgecolor','none'); alpha(Z,.05)
    for i = 1:size(aacenters,1)
        [V, Qout]  = shiftPoints(malloid2.V2P, axang2quat(aacenters(i,:)));
        Vh  = shiftPoints(malloid11.V2P, axang2quat(aacenters(i,:)));
        PlotMallow(malloid11,axang2quat(aacenters(i,:)),f,'red',.1,'none');
        atpi = find(abs(sqrt(sum(Vh.*Vh,2))-pi)<.0000001);
        scatter3(Vh(atpi,1),Vh(atpi,2),Vh(atpi,3),50,'r.');
        
        for tind = 1:size(malloid2.T2V,1)
            tri = malloid2.T2V(tind,:);
            
            aa1 = V(tri(1),:);
            aa2 = V(tri(2),:);
            aa3 = V(tri(3),:);
            
            q1 = axang2quat([aa1 norm(aa1)]);
            q2 = axang2quat([aa2 norm(aa2)]);
            q3 = axang2quat([aa3 norm(aa3)]);
            
            v1 = [q1(1) q2(1) q3(1)];
            %% search for a = [a1,a2,a3] s.t. [1 1 1]a' = 1, v1a' = 0, a > 0, a < 1
            % v1a' = 0 ensures the quat is a rotation by pi.
            nsp = null(v1);
            assert(sum(size(nsp)==[3,2])==2);
            dot2one = nsp'*[1;1;1];
            % find l,m s.t. dot2one'*[l;m]=1, 1 > nsp*[l;m] > 0
            % find specific solution by linprog with random objective.
            % Ax <= b Aeq x = beq
            Aeq = dot2one';
            beq = 1;
            A = [nsp;-nsp];
            b = [1;1;1;0;0;0];
            randobj = rand(1,2);
            [lm1,~,exitflag1] = linprog(randobj,A,b,Aeq,beq);
            [lm2,~,exitflag2] = linprog(-randobj,A,b,Aeq,beq);
            if exitflag1 == 1 && exitflag2 == 1 %-3 is unbounded
                % check solution is right
                assert(numel(lm1)~=0);
                BCcoords = nsp*lm1;
                assert(abs(v1*BCcoords)<.00000001);
                assert(abs(sum(BCcoords)-1)<.000000001);
                assert(sum(BCcoords>=-.000001)==3);
                assert(sum(BCcoords<=1.0000001)==3);
                assert(numel(lm2)~=0);
                BCcoords = nsp*lm2;
                assert(abs(v1*BCcoords)<.00000001);
                assert(abs(sum(BCcoords)-1)<.000000001);
                assert(sum(BCcoords>=-.000001)==3);
                assert(sum(BCcoords<=1.0000001)==3);
                
                QE1 = [q1' q2' q3']*nsp*lm1;
                QE2 = [q1' q2' q3']*nsp*lm2;
                
                if(norm(QE2-QE1)>.000001)
                    % compute solution line.
                    quatline = [linspace(QE1(1),QE2(1),lineres);linspace(QE1(2),QE2(2),lineres);linspace(QE1(3),QE2(3),lineres);linspace(QE1(4),QE2(4),lineres)]';
                    aa = quatToAxisAngle(quatline);
                    aa = aa(:,1:3).*aa(:,4);
                    plot3(aa(:,1),aa(:,2),aa(:,3),'g-','LineWidth',3);
                    
                    aaverts = [aaverts; aa];
                    solutionQEdges = [solutionQEdges; QE1' QE2'];
                end
            end
                
            
            
        end
        
        
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

