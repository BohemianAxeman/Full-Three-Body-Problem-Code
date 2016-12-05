% This function calculates the required spin rate and velocities of three bodies (in Barycenter
% Coordinates) such that they are tidally locked, given system angular momentum and relative
% distances between bodies, and their uniform mass and moment of inertia. System is assumed planar.

function [w_req, OutputBody] = wvHtidal(H, mass, I, r1, r2, r3)
for i=1:3
   OutputBody(i).v_req=zeros(1,3); OutputBody(i).v_vec=zeros(1,3); OutputBody(i).rB=zeros(1,3); 
end
RCM=(r1+r2+r3)/3; % Barycenter for system
r1B=r1-RCM; r2B=r2-RCM; r3B=r3-RCM; %Convert to Barycenter coordinates

w_req=H/(3*I+mass*(norm(r1B)^2+norm(r2B)^2+norm(r3B)^2)); %Required spin rate of bodies
v1=w_req*norm(r1B); v2=w_req*norm(r2B); v3=w_req*norm(r3B);  %Required tangential velocities
w_vec=[0,0,1]; %Rotation axis
OutputBody(1).v_req=v1; OutputBody(2).v_req=v2; OutputBody(3).v_req=v3; 

OutputBody(1).v_vec=v1*cross(w_vec,r1B)/norm(cross(w_vec,r1B)); %Direction of veloicty vector for body 1
OutputBody(2).v_vec=v2*cross(w_vec,r2B)/norm(cross(w_vec,r2B)); %Direction of veloicty vector for body 2
OutputBody(3).v_vec=v3*cross(w_vec,r3B)/norm(cross(w_vec,r3B)); %Direction of veloicty vector for body 3

OutputBody(1).rB=r1B; OutputBody(2).rB=r2B; OutputBody(3).rB=r3B; 

end
