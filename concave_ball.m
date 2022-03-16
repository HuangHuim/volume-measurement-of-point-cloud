close all;clear;clc;
%% generate the point cloud
r=1; % the ridus of the sphere
a=0;b=3;c=0;d=4;
l1=abs(b-a);l2=abs(d-c);
n=10; % partition number
% Equidistant partition
x=linspace(a,b,l1*n+1);
y=linspace(c,d,l2*n+1);
[X,Y]=meshgrid(x,y);
[M,N]=size(X);
Z=zeros(M,N);
a1=(-2)*sqrt(2);
for i=1:M  % y label
    for j=1:N
        t1=2*X(i,j);
        if Y(i,j)<1/2
            if  Y(i,j)>t1
                Z(i,j)=2*a1*X(i,j);
            elseif Y(i,j)<-1*t1+6
                Z(i,j)=a1*Y(i,j);
            else
                Z(i,j)=(-2)*a1*X(i,j)+6*a1;
            end
        elseif Y(i,j)<7/2
            if X(i,j)<1/4
                Z(i,j)=2*a1*X(i,j);
            elseif X(i,j)<11/4
                t2=(X(i,j)-3/2)^2+(Y(i,j)-2)^2;
                if t2<=1
                    Z(i,j)=sqrt(1-t2)-sqrt(2);
                else
                    Z(i,j)=-sqrt(2);
                end
            else
                Z(i,j)=(-2)*a1*X(i,j)+6*a1;
            end
        else
            if Y(i,j)<-1*t1+4
                Z(i,j)=2*a1*X(i,j);
            elseif Y(i,j)>t1-2
                Z(i,j)=2*sqrt(2)*Y(i,j)-8*sqrt(2);
            else
                Z(i,j)=(-2)*a1*X(i,j)+6*a1;
            end
        end
    end
end


% Z=-1.*sqrt(1-X.^2-Y.^2);
Z=real(Z);
figure;surf(X,Y,Z,'EdgeColor','none');
xlabel('X');ylabel('Y');zlabel('Z');
% ball_01=[X,Y,Z];
%% the second-order estimation
delt_r=1/n;
s=delt_r*delt_r;
V=zeros(M,N);
for i=1:M-1
    for j=1:N-1
        f1=Z(i,j);
        f2=Z(i+1,j);
        f3=Z(i,j+1);
        f4=Z(i+1,j+1);
        average_h=(f1+f2+f3+f4)/4;
        V(i,j)=s*average_h;
    end
end
Volume=sum(sum(V));

% the second-order estimation
% filter 6*6
H=[0 0 4 4 0 0;0 3 -7 -7 3 0;4 -7 3 3 -7 4;4 -7 3 3 -7 4;0 3 -7 -7 3 0;0 0 4 4 0 0];
figure;mesh(H);
% mf=6; % the size of the filter
% half_mf=5*mf/2;
Hf=filter2(H,Z,'same');
% M_low=half_mf;  M_high=M-half_mf-1;
% N_low=half_mf;  N_high=N-half_mf-1;
% h2V=Hf(M_low:M_high,N_low:N_high);
h2V=Hf(1:M-1,1:N-1);

dividend=4*12*4;  %four direction, four points, 12 is the common dividend
V_Hf=sum(sum(Hf))*s/dividend;
V222=sum(sum(h2V))*s/dividend;

Hf_full=filter2(H,Z,'full');
V_full=sum(sum(Hf_full))*s/dividend;

figure;surf(Hf,'EdgeColor','none');