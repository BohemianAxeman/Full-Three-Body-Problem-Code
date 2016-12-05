
function BodyOut = calcForces2(CurrState, G, N, h, breakOut)

global MU DU TU
BodyOut=CurrState; %used to transfer static body properties

% Organize body properties into easy to work with vector arrays
    r_vec(1:N,1:3)=CurrState.r(1:N,1:3);
    v_vec(1:N,1:3)=CurrState.v(1:N,1:3);
    psi_vec(1:N,1)=CurrState.psi(1:N,1);
    w_vec(1:N,1:3)=CurrState.w(1:N,1:3);
    contact(1:N,1:N)=CurrState.contact(1:N,1:N);
    m(1:N)=CurrState.m(1:N);
    R(1:N)=CurrState.R(1:N);
    Ivec(1:N,1:3)=CurrState.Ivec(1:N,1:3);
    
totForce = zeros(N,3);
% totForce2 = zeros(N,3);

% %%%%%%%%%%%%%%%%%%%% Calculate Gravitational Forces %%%%%%%%%%%%%%%%%%%%
for k = 1:N-1
    for j = k+1:N
       rkj = [r_vec(j,1) - r_vec(k,1), r_vec(j,2) - r_vec(k,2), r_vec(j,3) - r_vec(k,3)];
       Frel = m(k)*m(j)*(rkj)./norm(rkj)^3; % Force felt by mass k due to j
       totForce(k,1:3) = totForce(k,1:3) + Frel; 
       totForce(j,1:3) = totForce(j,1:3) - Frel; 
    end
end
totForce = totForce*G;


%%%%%%%%%%%%%%%%%%%% Calculate Contact Spring Forces %%%%%%%%%%%%%%%%%%%%
if any(any(contact))
    ContactPoints = breakOut.ContactPoints;
    ke = 1e13 * TU^2/MU; %Spring elastic constant, [kN/km]
    kd = 1e13 * TU/MU; %Spring damping constant, [kN.s/km]
for nc = 1:breakOut.nc
    if ~breakOut.BC(nc) %Bodies not separating
        k=ContactPoints(nc).A;
        j=ContactPoints(nc).B;
        rkj = r_vec(j,1:3) - r_vec(k,1:3);
        n = rkj/norm(rkj);
        vkp = v_vec(k,1:3);% + cross(w_vec(k,1:3), R(k)*n); %Velocity of contact point on k
        vjp = v_vec(j,1:3);% + cross(w_vec(j,1:3), -R(k)*n); %Velocity of contact point on j
        vkj = vjp - vkp; %Relative contact point velocity
        z = R(k)+R(j) - norm(rkj); %Penetration distance
        zdot = dot(vkj, n);
        Fe = ke*z*n; %Elastic spring force
        Fd = -kd*zdot*n; %Damper spring force

        totForce(j,1:3) = totForce(j,1:3) + Fe + Fd;
        totForce(k,1:3) = totForce(k,1:3) - Fe - Fd;

    end
end
end

%%%%%%%%%%%%%%%%%%%% Calculate Ellipsoidal Gravity & Tidal Torques %%%%%%%%%%%%%%%%%%%%
psi_ddot_sum = zeros(1,N); %Initialize angular acceleration scalar
if true
    epsilon = 5e-3; % Specific tidal energy dissipation constant
    klove = 1e-5; %.0042; % Tidal love number (.03 = moon)
    rho = 2e12 * (DU^3/MU); %[kg/km^3], density of asteroid
for k=1:N
%     ddot_Xsum=0; ddot_Ysum=0;
    for j=1:N
        if j==k; continue; end
        rkj = r_vec(j,1:3) - r_vec(k,1:3);
        n = rkj/norm(rkj);
        theta = atan2(rkj(2),rkj(1)); %Angle between connecting body axis and X direction
        % phi = psi - theta
            phik = psi_vec(k)-theta; %Angle from rjk to body buldge k
            phij = psi_vec(j)-theta; %Angle from rjk to body buldge j
        % Calculate phi_dot & theta_dot
            theta_dot=((r_vec(j,1)-r_vec(k,1))*(v_vec(j,2)-v_vec(k,2))-((r_vec(j,2)-r_vec(k,2))*...
                (v_vec(j,1)-v_vec(k,1))))/((r_vec(j,1)-r_vec(k,1))^2+(r_vec(j,2)-r_vec(k,2))^2);
            phi_dot=w_vec(k,3)-theta_dot;
%             if norm(phi_dot) <1e-9
%                 phi_dot = 0;
%             end
        % Calculate the tidal torque, Gamma
            Gamma = sign(phi_dot)*3/2*klove*(3/(4*pi*rho))^2*G*m(k)^2*m(j)^2/...
                (norm(rkj)^6*R(k))*sin(2*epsilon); %Tidal torque
        % Must check if phi_dot is passing near 0
            delta = norm(Gamma/Ivec(k,3)*sqrt(6/(pi*G*rho))); %Small Characteristic spin rate
            if abs(phi_dot) < delta %Must use modified tidal torque equation
                Gamma=sqrt(pi*G*rho/6)*Ivec(k,3)*phi_dot;
%                 if norm(Gamma)<1e-6 %Gamma probably supposed to be 0, account for numerical errors
%                     Gamma = 0;
%                 end
            end
            
            
        psi_ddot = 1/Ivec(k,3)*(-3*G/2*(Ivec(k,2)-Ivec(k,1))*m(j)*sin(2*phik)/norm(rkj)^3 - Gamma);
        psi_ddot_sum(k) = psi_ddot_sum(k) + psi_ddot/2;
           
        % Clean up small numerical inaccuracies
        %         if norm(psi_ddot_sum(k))<eps; 
        %             psi_ddot_sum(k)=0; 
        %         end
        
        %Calculate Corresponding change in orbital energy
        F_tidal = Gamma/(norm(rkj)); %Effective tangential force due to tidal torque on body j
%         if phi_dot > 0 && sum(contact(k,:))>0 
%             %Don't want bodies to separate if in contact
%             F_tidal = F_tidal; %*.947;
%         end
        n_tan = cross([0,0,1],n); %theta_dot direction
        totForce(j,:) = totForce(j,:) + F_tidal*n_tan; %Apply force in theta_dot direction
        
%         % **** Sum Effects of gravity potential on kth body due to jth body ****
%             pr_px=(r_vec(k,1)-r_vec(j,1)) / norm(rkj); %Partial of rjk wrt x_k
%             pr_py=(r_vec(k,2)-r_vec(j,2)) / norm(rkj); %Partial of rjk wrt y_k
%             ptheta_px=(r_vec(j,2)-r_vec(k,2)) / norm(rkj)^2; pphi_px=-ptheta_px; % Partial phi wrt x_n
%             ptheta_py=(r_vec(k,1)-r_vec(j,1)) / norm(rkj)^2; pphi_py=-ptheta_py; % Partial phi wrt y_n
%             b=sum(Ivec(k,1:3)) + sum(Ivec(j,1:3)); %Sum of kth and jth bodies moments of intertia
%             c=Ivec(k,1)+Ivec(k,2)+Ivec(j,1)+Ivec(j,2);
%             d=Ivec(k,2)-Ivec(k,1); %This should be the same for j body, else need to change code partially
%             % ddot_x = -1/m * P(V)/P(x)   ddot_y = -1/m * P(V)/P(y)
%             ddot_X=-1/norm(rkj)^2 * (pr_px* (m(k)+ 3/(2*norm(rkj)^2)* (b-3/2*( c-d*cos(phik)-d*cos(phij))))...
%                 + 3*d/(4*norm(rkj))*pphi_px * (sin(phik)+sin(phij)));
%             ddot_Y=-1/norm(rkj)^2 * (pr_py* (m(k)+ 3/(2*norm(rkj)^2)* (b-3/2*( c-d*cos(phik)-d*cos(phij))))...
%                 + 3*d/(4*norm(rkj))*pphi_py * (sin(phik)+sin(phij)));
% 
%             ddot_Xsum=ddot_Xsum + ddot_X;
%             ddot_Ysum=ddot_Ysum + ddot_Y;

    end %End of j for loop
    
         
%     totForce2(k,1:3) = G*[ddot_Xsum, ddot_Ysum, 0];
end %End of k for loop

end %Turn tidal torques on/off


%%%%%%%%%%%%%%%%%%%% Calculate Geometric Contact Forces %%%%%%%%%%%%%%%%%%%%
% extF = totForce;
% for i=1:N
%    extT(i,1:3)=[0,0 Ivec(i,3)*psi_ddot_sum(i)]; 
% end
% 
% 
% if any(any(contact))
%     augcontact = contact;
%     for i=1:breakOut.nc
%         if breakOut.BC(i) %Broke contact, don't apply contact forces
%             A = breakOut.ContactPoints(i).A; B = breakOut.ContactPoints(i).B;
%             augcontact(A,B) = 0; augcontact(B,A) = 0;
%         end
%     end
%     augState=CurrState; augState.contact = augcontact;
%     %Bodies are in contact, must calculate contact forces
%     [add_F, x, exitFlag, output, returnPacket]=contactForces(augState, extF, extT, N);
%     if exitFlag~=1
%        disp('Contact force calculation failed') 
%        disp('Iteration #: ') 
%        %disp(output.iterations)
%        disp('Final vector value reached: ')
%        disp(x)
%        disp('Time of failure: ')
%        disp(currTime)
%        while true
%         switch input('Select an action: (0 to quit, 1 to continue, 2 for error readout, 3 for debug) ')
%            case 0
%                error('User exited due to failure of contact force convergence')
%            case 1 
%                disp('Continued'), break;
%            case 2
%                disp(output.message)
%            case 3
%                disp('Entered debug mode')
%                dbstop if warning
%                warning('Debug mode selected')
%         end
%        end
%     end
%     
%     %Add in contact forces to net acceleration vector
%     for i=1:N
%        totForce(i,:)=totForce(i,1:3) + add_F(i,1:3);
% %        for k=1:3
% %           if norm(BodyOut.v(i,k))<eps
% %               BodyOut.v(i,k)=0;
% %           end
% %        end
%     end
% end

%%%%%%%%%%%%%%%%%%%% Return net acceleration %%%%%%%%%%%%%%%%%%%%
for k=1:N
   BodyOut.r(k,1:3) = CurrState.v(k,1:3);
   BodyOut.v(k,1:3) = totForce(k,1:3)/m(k);
   BodyOut.psi(k,1) = CurrState.w(k,3);
   BodyOut.w(k,1:3) = [0,0,psi_ddot_sum(k)];
end


end