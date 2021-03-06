function X = Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1, x2, x3, X0)
%% Nonlinear_Triangulation
% Refining the poses of the cameras to get a better estimate of the points
% 3D position
% Inputs: 
%     K - size (3 x 3) camera calibration (intrinsics) matrix
%     x
% Outputs: 
%     X - size (N x 3) matrix of refined point 3D locations
N=size(X0,1);
X=size(N,3);
X=X0;
for i=1:10
    for j=1:N
      X(j,:)=Single_Point_Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1(j,:), x2(j,:), x3(j,:), X(j,:));
    end
end
end

function X = Single_Point_Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1, x2, x3, X0)
J=[Jacobian_Triangulation(C1, R1, K, X0) ...
    Jacobian_Triangulation(C2, R2, K, X0) ...
    Jacobian_Triangulation(C3, R3, K, X0)]';
b=[x1 x2 x3]'; 
uvw1=K*R1*[X0'-C1];
uvw2=K*R2*[X0'-C2];
uvw3=K*R3*[X0'-C3];
F=[uvw1(1)/uvw1(3); uvw1(2)/uvw1(3) ;uvw2(1)/uvw2(3); ...
    uvw2(2)/uvw2(3) ;uvw3(1)/uvw3(3); uvw3(2)/uvw3(3)];
del_x=pinv(J'*J)*J'*(b-F);
X=del_x'+X0;


end

function J = Jacobian_Triangulation(C, R, K, X)
uvw=K*R*[X'-C];
f=K(1,1);px=K(1,end);py=K(2,end);
deluBydeldx=[f*R(1,1)+px*R(3,1) f*R(1,2)+px*R(3,2)  f*R(1,3)+px*R(3,3)];
delvBydeldx=[f*R(2,1)+py*R(3,1) f*R(2,2)+py*R(3,2)  f*R(2,3)+py*R(3,3)];
delwBydeldx=[R(3,1) R(3,2) R(3,3)];
delfBydeldx=[((uvw(3)*deluBydeldx)-(uvw(1)*delwBydeldx))/(uvw(3)^2);...
    ((uvw(3)*delvBydeldx)-(uvw(2)*delwBydeldx))/(uvw(3)^2)];
J=delfBydeldx';
end
