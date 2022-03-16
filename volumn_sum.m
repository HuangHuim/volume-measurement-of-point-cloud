% load CoorP.mat;
close all;
clear;clc
load SinglePoint.mat;
load ProCoor.mat;    %load the project point in the road plain
load reserved_Tri.mat;
q=SingleHole;
p=ProCoor;
T3=reserved_Tri;    % simiplify mark
nn=size(reserved_Tri,1);
V_sum=zeros(nn,1);
for j=1:nn
    V1=abs(dot(q(T3(j,1),:)-p(T3(j,1),:),cross(p(T3(j,2),:)-p(T3(j,1),:),p(T3(j,3),:)-p(T3(j,1),:))));
    V2=abs(dot(q(T3(j,2),:)-p(T3(j,2),:),cross(p(T3(j,1),:)-p(T3(j,2),:),p(T3(j,3),:)-p(T3(j,2),:))));
    V3=abs(dot(q(T3(j,3),:)-p(T3(j,3),:),cross(p(T3(j,2),:)-p(T3(j,3),:),p(T3(j,1),:)-p(T3(j,3),:))));
    V_sum(j)=(V1+V2+V3)/6;
end
V=sum(V_sum);
    

S_sum=zeros(nn,1);
for j=1:nn
    S1=norm(cross(p(T3(j,2),:)-p(T3(j,1),:),p(T3(j,3),:)-p(T3(j,1),:)));
    S2=norm(cross(p(T3(j,1),:)-p(T3(j,2),:),p(T3(j,3),:)-p(T3(j,2),:)));
    S3=norm(cross(p(T3(j,2),:)-p(T3(j,3),:),p(T3(j,1),:)-p(T3(j,3),:)));
    S_sum(j)=(S1+S2+S3)/6;
end
S=sum(S_sum);

width=max(ProCoor(:,1))-min(ProCoor(:,1));
length=max(ProCoor(:,2))-min(ProCoor(:,2));
hight=max(SingleHole(:,3))-min(SingleHole(:,3));
