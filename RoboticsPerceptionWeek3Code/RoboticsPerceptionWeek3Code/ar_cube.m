function [proj_points, t, R] = ar_cube(H,render_points,K)
%% ar_cube
% Estimate your position and orientation with respect to a set of 4 points on the ground
% Inputs:
%    H - the computed homography from the corners in the image
%    render_points - size (N x 3) matrix of world points to project
%    K - size (3 x 3) calibration matrix for the camera
% Outputs: 
%    proj_points - size (N x 2) matrix of the projected points in pixel
%      coordinates
%    t - size (3 x 1) vector of the translation of the transformation
%    R - size (3 x 3) matrix of the rotation of the transformation
% Written by Stephen Phillips for the Coursera Robotics:Perception course

% YOUR CODE HERE: Extract the pose from the homography

% YOUR CODE HERE: Project the points using the pose
t=H(:,3)/norm(H(:,1));
if t(3)<0
    H=-H;
end
t=H(:,3)/norm(H(:,1));
R=[H(:,1) H(:,2) cross(H(:,1),H(:,2))];
[U,S,Vt]=svd(R);
R=U*[1 0 0;0 1 0;0 0 det(U*Vt)]*Vt;
render_points(:,end+1)=ones(size(render_points,1),1);
proj_points=K*[R t]*render_points';
proj_points=proj_points./proj_points(3,:);
proj_points=proj_points(1:2,:);
proj_points=proj_points';


end
