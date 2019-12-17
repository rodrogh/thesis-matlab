%Module MECH3460
%Transformation Matrices

%Constants
clear all;
clc;
% alpha = sym([pi/2 0 -pi/2 pi/2 -pi/2 0]);
% a = sym([0 0.432 0.20 0 0 0]);
% theta = sym('theta',[1 6]);
% d = sym([0 0 0.150 0.432 0 0]);
alpha = sym([ pi/2 0 0 ]);
theta = sym([ sym('theta1') sym('theta2') sym('theta3')]);
% theta = sym([pi pi/2 sym('theta3')]);
a = sym([ sym('L1') sym('L2') 0]);
% a = sym([ sym('L1') sym('L2')*sin(theta(3)) 0 ]);
d = sym([0 0 0]);
temp=1;
for i=1:3
    T_sym(:,:,i) = [cos(theta(i))  -cos(alpha(i))*sin(theta(i))   sin(alpha(i))*sin(theta(i))   a(i)*cos(theta(i));...
                sin(theta(i))  cos(alpha(i))*cos(theta(i))   -sin(alpha(i))*cos(theta(i))   a(i)*sin(theta(i));...
                0               sin(alpha(i))                   cos(alpha(i))                  d(i);...
                0               0                                0                               1];
    temp = temp*T_sym(:,:,i);
end
disp(T_sym);
disp(temp);
% alpha = [90 0 -90 90 -90 0];
% a = [0 0.432 0.20 0 0 0];
% theta = [0 0 0 0 -30 0];
% d = [0 0 0.150 0.432 0 0];
% temp = 1;
% for i=1:6
%     T(:,:,i) = [cosd(theta(i))  -cosd(alpha(i))*sind(theta(i))   sind(alpha(i))*sind(theta(i))   a(i)*cosd(theta(i));...
%                 sind(theta(i))  cosd(alpha(i))*cosd(theta(i))   -sind(alpha(i))*cosd(theta(i))   a(i)*sind(theta(i));...
%                 0               sind(alpha(i))                   cosd(alpha(i))                  d(i);...
%                 0               0                                0                               1];
%     temp = temp*T(:,:,i);   
% end
% T0_6 = temp;
% disp(T);
% disp(T0_6);
