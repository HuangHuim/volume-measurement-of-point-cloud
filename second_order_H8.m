function [V_vector,E_vector]=second_order_H8(Z,s,SS)

H8=[0 0 2 2 2 2 0 0;0 3 -2 -5 -5 -2 3 0;2 -2 -4 0 0 -4 -2 2;2 -5 0 7 7 0 -5 2;...
    2 -5 0 7 7 0 -5 2;2 -2 -4 0 0 -4 -2 2;0 3 -2 -5 -5 -2 3 0;0 0 2 2 2 2 0 0];
% figure;surf(H8);
Hf_8=filter2(H8,Z,'same');
averageH_8=sum(sum(Hf_8))/(M*N);
V_Hf_averageH_8=averageH_8*SS/dividend;

figure;subplot(1,2,1);surf(X,Y,Hf_8,'EdgeColor','none');
axis tight;axis off;
% axis([-1,1,-1,1]);
subplot(1,2,2);imagesc([-1;1],[-1;1],Hf_8);
axis tight;axis off;