function test
close all;

Amax = pi+2;
MaxIter = 100000;
n = 300;
recton = true;
h = 1;

L = sparse(n,n);
L(1:n+1:end)=1;
L(2:n+1:end)=-1/2;
L(n+1:n+1:end)=-1/2;
L(1,n)=-1/2; L(n,1)=-1/2;
L = kron(L,speye(2));
P = L;

A = sparse(2*n,2*n);
A(sub2ind(size(A), 1:2:2*(n-1), 4:2:2*n)) = 1;
A(end-1,2) = 1;
A(sub2ind(size(A), 3:2:2*n, 2:2:2*(n-1))) = -1;
A(1,end) = -1;
A = -A*.5;
A = (A + A') /2;
eigs(A)

x0 = [1:2*(n+1)]';
x0(1:2:end) = .999*cos(linspace(0,2*pi,n+1));
x0(2:2:end) = .999*sin(linspace(0,2*pi,n+1));
x0 = x0(1:2*n);

% ellipse
R1 = 1;
R2 = Amax/(pi*R1);
x0 = [1:2*(n+1)]';
x0(1:2:end) = R2*cos(linspace(0,2*pi,n+1));
x0(2:2:end) = R1*sin(linspace(0,2*pi,n+1));
x0 = x0(1:2*n);
assert(abs(abs(x0'*A*x0)-Amax)<.001);
eperim = 2*pi*sqrt((R1^2+R2^2)/2);
assert(abs(perimeter(x0)-eperim)<1);

%x0 = 100*(rand(2*n,1)-.5); % xyxyxyxyxyxyxy....xyxyxyxy

% x0 = [0,0,0,1,1,1,1,0]';

Perim = x0'*P*x0
Area = x0'*A*x0

assert(abs(abs(Area) - polyarea(x0(1:2:end),x0(2:2:end)))<.00001);
assert(abs(Perim*2 - norm(x0 - circshift(x0,2))^2)<.000001);

%% start optimization

fun = @(x)quadobj(x,P);
constr = @(x)quadconstr(x,-A,Amax);
lb = [1:2*n]; lb(2:2:end) = -1; lb(1:2:end) = -10000;
ub = [1:2*n]; ub(2:2:end) = 1; ub(1:2:end) = 10000;

find(x0>=ub')
find(x0<=lb')

A0 = x0'*A*x0
P0 = x0'*P*x0

options = optimoptions(@fmincon,'Algorithm','interior-point',...
    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
    'display','iter','CheckGradients',false,'FiniteDifferenceType','central',...
    'SubproblemAlgorithm','cg','HessianMultiplyFcn',@(x,l,v)HessMultFcn(x,l,v,fun,constr),'MaxPCGIter',50);
options.MaxFunctionEvaluations = MaxIter;
options.MaxIterations = MaxIter;



if(recton)
    [x,fval,eflag,output,lambda] = fmincon(fun,x0,[],[],[],[],lb',ub',constr,options);
else
    [x,fval,eflag,output,lambda] = fmincon(fun,x0,[],[],[],[],[],[],constr,options);
end

objinit = fun(x0)
[~, constrinit] = constr(x0)
Pinit = perimeter(x0)

objfinal = fun(x)
[~, constrfinal] = constr(x)
Pfinal = perimeter(x)

Area = x'*A*x
Perim = x'*P*x

V0 = reshape(x0,2,[])';
V = reshape(x,2,[])';
hold on; axis equal
plot([-1 1],[1 1],'g-');
plot([-1 1],[-1 -1],'g-');

plot(V(:,1),V(:,2),'k-');
plot(V(:,1),V(:,2),'r.');

plot(V0(:,1),V0(:,2),'g-');
plot(V0(:,1),V0(:,2),'g.');

end


function W = HessMultFcn(x,lambda,v,obj,constr)
    
    [a,b] = obj(v);
    [c,d,e,f]=constr(v);
    W = f*lambda.eqnonlin + b;
end

function [y,grady] = quadobj(x,Q)
    y = x'*Q*x;
    
    if nargout > 1
        grady = 2*Q*x;
    end
end

function [y,yeq,grady,gradyeq] = quadconstr(x,Q,k)
    y = []; grady = [];
    yeq = x'*Q*x - k;
    
    if nargout > 2
        gradyeq = 2*Q*x;
    end
end


