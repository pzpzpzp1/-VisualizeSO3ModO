figure; hold on; axis equal; n = 100000;
scatter3(0,0,0,'r');
aa = [0,0,pi/2];
scatter3(aa(1),aa(2),aa(3),'r');
z1 = Bisector(aa,n);
scatter3(z1(:,1),z1(:,2),z1(:,3),1,'r');

aa = (2*pi/3)*[1,1,1]/sqrt(3);
scatter3(aa(1),aa(2),aa(3),'r');
z2 = Bisector(aa,n);
scatter3(z2(:,1),z2(:,2),z2(:,3),1,'g');

aa = (2*pi/3)*[-1,1,1]/sqrt(3);
scatter3(aa(1),aa(2),aa(3),'r');
z2 = Bisector(aa,n);
scatter3(z2(:,1),z2(:,2),z2(:,3),1,'g');

aa = (2*pi/3)*[1,-1,1]/sqrt(3);
scatter3(aa(1),aa(2),aa(3),'r');
z2 = Bisector(aa,n);
scatter3(z2(:,1),z2(:,2),z2(:,3),1,'g');

aa = (2*pi/3)*[-1,-1,1]/sqrt(3);
scatter3(aa(1),aa(2),aa(3),'r');
z2 = Bisector(aa,n);
scatter3(z2(:,1),z2(:,2),z2(:,3),1,'g');
