clear;
n = 100;
p = 3;
theta = linspace(0,pi,n);
phi = linspace(0,2*pi,n);

for i=1:n
    for j=1:n
        X(i,j) = sin(theta(i))*cos(phi(j));
        Y(i,j) = sin(theta(i))*sin(phi(j));
        Z(i,j) = cos(theta(i));
    end
end

v = {X.*X, X.*Y, X.*Z, Y.*Y, Y.*Z, Z.*Z};
r = zeros(size(X));

M = randn(6,6); % random semidefinite matrix
M = M'*M;

for i=1:6
    for j=1:6
        r = r + M(i,j)*v{i}.*v{j};
    end
end

subplot(1,2,1);
surf(real(X).*r,real(Y).*r,real(Z).*r,real(r),'edgecolor','none');
axis equal;
cameratoolbar;
title('Random quartic SOS');

subplot(1,2,2);
r = max(r(:)) - r;
surf(real(X).*r,real(Y).*r,real(Z).*r,real(r),'edgecolor','none');
axis equal;
cameratoolbar;
title('Flipped max and min');