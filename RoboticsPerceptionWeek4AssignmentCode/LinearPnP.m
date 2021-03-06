function [C, R] = LinearPnP(X, x, K)
%% LinearPnP
% Getting pose from 2D-3D correspondences
% Inputs:
%     X - size (N x 3) matrix of 3D points
%     x - size (N x 2) matrix of 2D points whose rows correspond with X
%     K - size (3 x 3) camera calibration (intrinsics) matrix
% Outputs:
%     C - size (3 x 1) pose transation
%     R - size (3 x 1) pose rotation
%
% IMPORTANT NOTE: While theoretically you can use the x directly when solving
% for the P = [R t] matrix then use the K matrix to correct the error, this is
% more numeically unstable, and thus it is better to calibrate the x values
% before the computation of P then extract R and t directly
N=size(x,1);
P=zeros(3,4);
x(:,end+1)=ones(N,1);
x=x';
X(:,end+1)=ones(N,1);
xc=pinv(K)*x;
t=zeros(3,1);
for i=1:N
    A(3*i-2,:)=[0 0 0 0 -X(i,:) xc(2,i)*X(i,:)];
     A(3*i-1,:)=[X(i,:) 0 0 0 0 -xc(1,i)*X(i,:)];
     A(3*i,:)=[-xc(2,i)*X(i,:) xc(1,i)*X(i,:) 0 0 0 0];
end
[U,D,Vt]=svd(A);
P=(reshape(Vt(:,end),4,3))';
R=P(:,1:3);
[U,D,Vt]=svd(R);
if(det(U*Vt')>0)
    R=U*Vt';
    t=(P(:,4))/D(1,1);
else
    R=-U*Vt';
    t=-(P(:,4))/D(1,1);
end
C=-R'*t;






