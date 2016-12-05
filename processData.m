
%%%%%%%%%%%%%%%% State Analysis and Data Visualization %%%%%%%%%%%%%%%% 
global TU MU DU
G=6.67384e-20; % km^3/(kg-s^2)
N = 3;

% Post processing of the output from mainSimulation
close all

%%%%% Need to elminate trailing 0 elements of time and body vectors
if tvec(end) == 0
    indDataEnd = find(tvec(2:end)==0,1) + 1; %add 1 to index beacuse we leave out first element of tvec (tvec(1)=0)
else
    indDataEnd = length(tvec) + 1;
end
    t1 = tvec(1:indDataEnd-1)*TU;

for k=1:N
    rplot(k).r(:,1:3) = BodyOut(k).r(1:indDataEnd-1,1:3);
    rplot(k).v(:,1:3) = BodyOut(k).v(1:indDataEnd-1,1:3);
    rplot(k).w(:,1:3) = BodyOut(k).w(1:indDataEnd-1,1:3);
    rplot(k).psi(:,1) = BodyOut(k).psi(1:indDataEnd-1,1);
    rplot(k).mass = BodyOut.mass; rplot(k).I=BodyOut.I; rplot(k).R=BodyOut.R;
    rplot(k).Ivec(1:3)=BodyOut(k).Ivec(1:3); 
    rplot(k).contact(:,1:N)=BodyOut(k).contact(1:indDataEnd-1,:);
end
%%%%%%

radius=rplot(1).R*DU;

%Preallocate Arrays
arraySize=length(t1);
r=zeros(arraySize, 3, N); v=zeros(arraySize, 3, N); w=zeros(arraySize, 3, N); psi=zeros(arraySize, 1, N);

%Pull out useful information for further computing
for k=1:N
   r(:,1:3,k) = (DU) * rplot(k).r(:,1:3); %position vector: row (time step), column (x,y,z), page (body #)
   v(:,1:3,k) = (DU/TU) * rplot(k).v(:,1:3); %velocity vector: row (time step), column (x,y,z), page (body #)
   w(:,1:3,k) = (1/TU) * rplot(k).w(:,1:3); %rad/s
   psi(:,1,k) = rplot(k).psi(:,1); %rad
   mass(k) = (MU) * rplot(k).mass; I(k) = (MU*DU^2) * rplot(k).I; R(k) = (DU) * rplot(k).R;
   Ivec(k,1:3) = (MU*DU^2) * rplot(k).Ivec(1:3); Isum(k) = (MU*DU^2) * sum(Ivec(k,1:3));
   contacts(:,1:N,k) = rplot(k).contact(:,1:N);
end

T=zeros(length(t1),1); U=zeros(length(t1),1); E=zeros(length(t1),1); H=zeros(length(t1),3);
Hmag=zeros(length(t1),1); Ip=zeros(length(t1),1); Ih=zeros(length(t1),1); phi=zeros(length(t1),N,N);
Rcm=zeros(length(t1),3); RcmMAG=zeros(length(t1),1); Ttrans=zeros(length(t1),1); Trot=zeros(length(t1),1);
rMAG=zeros(length(t1),2); vMAG=zeros(length(t1),2); Hrv=zeros(length(t1),3); HwI=zeros(length(t1),3);
Trel=zeros(length(t1),1); Hrel=zeros(length(t1),3); Hrvmag=zeros(length(t1),1); HwImag=zeros(length(t1),1);
theta=zeros(length(t1),N,N); phi_dot=zeros(length(t1),N,N); theta_dot=zeros(length(t1),N,N);
B1RM=@(phi) [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1]; vCM=zeros(length(t1),3);
Trel_CM=zeros(length(t1),1); Hrel_CM=zeros(length(t1),3); vrel_CM=zeros(length(t1),3,N);
vpar=zeros(length(t1),N); vperp=zeros(length(t1),N); Ediff=zeros(length(t1)-1,1); Hdiff=zeros(length(t1)-1,1);
Hrv_multi=zeros(length(t1),3,N);
for j=1:length(t1)
    for k=1:N
        %Compute Center of Mass at each time step
        Rcm(j,1:3)=Rcm(j,1:3)+(mass(k)*r(j,1:3,k))/sum(mass);
    end
    if j==1
        continue
    end
    vCM(j,1:3)=(Rcm(j,:)-Rcm(j-1,:))/(t1(j)-t1(j-1));
end
for j=1:length(t1)
    for k=1:N
        %Compute and sum kinetic energy and angular momentum for each body
        Ttrans(j)=Ttrans(j) + 1/2*mass(k)*dot(v(j,1:3,k),v(j,1:3,k));
        Trot(j)=Trot(j) + 1/2*w(j,1:3,k)*diag(Ivec(k,1:3))*w(j,1:3,k)';
        Hrv(j,1:3)=Hrv(j,1:3) + mass(k)*cross(r(j,1:3,k),v(j,1:3,k));
        HwI(j,1:3)=HwI(j,1:3) + w(j,1:3,k)*diag(Ivec(k,1:3));
        Hrv_multi(j,1:3,k)=Hrv_multi(j,1:3,k) + mass(k)*cross(r(j,1:3,k),v(j,1:3,k));
        Ip(j)=Ip(j)+mass(k)*dot(r(j,1:3,k),r(j,1:3,k));

        %Compute r and v magnitudes
        rMAG(j,k)=norm(r(j,1:3,k)); %km,   row (time step), column (body #)
        vMAG(j,k)=norm(v(j,1:3,k)); %km/s, row (time step), column (body #)
        
        %Compute T and H using Rcm to update barycenter
        vrel_CM(j,1:3,k)=v(j,1:3,k)-vCM(j,1:3);
        Trel_CM(j)=Trel_CM(j) + 1/2*mass(k)*dot(vrel_CM(j,1:3,k),vrel_CM(j,1:3,k)) + 1/2*w(j,1:3,k)*diag(Ivec(k,1:3))*w(j,1:3,k)';
        Hrel_CM(j,1:3)=Hrel_CM(j,1:3) + mass(k)*cross(r(j,1:3,k)-Rcm(j,:),vrel_CM(j,1:3,k)) + w(j,1:3,k)*diag(Ivec(k,1:3));
        
%         %Compute rotating Body1 coordinates of each Body
%         rB1(j,:,k)=B1RM(psi(j,3,1))*(r(j,:,k)-r(j,:,1))';
    end
        RcmMAG(j)=norm(Rcm(j,1:3));
        T(j)=Ttrans(j)+Trot(j);
        H(j,1:3)=Hrv(j,1:3)+HwI(j,1:3);
    
    % Compute relative body spin angle and spin rate, phi & phi_dot
    for k=1:N 
        for m=1:N
            if k==m; continue; end
            theta(j,k,m)=atan((r(j,2,m)-r(j,2,k)) / (r(j,1,m)-r(j,1,k))); theta(j,m,k)=theta(j,k,m);
            phi(j,k,m)=psi(j,1,k) - theta(j,k,m); phi(j,m,k)=psi(j,1,m) - theta(j,k,m);
            theta_dot(j,k,m)=((r(j,1,m)-r(j,1,k))*(v(j,2,m)-v(j,2,k))-((r(j,2,m)-r(j,2,k))*(v(j,1,m)-v(j,1,k))))/...
            ((r(j,1,m)-r(j,1,k))^2+(r(j,2,m)-r(j,2,k))^2);
            phi_dot(j,k,m)=w(j,3,k)-theta_dot(j,k,m);
        end
    end
    % Compute potential energy of system, avoiding any self potentials
    Uint=zeros(N,N);
    for k=1:N-1
        for m=k+1:N
            rrel=norm(r(j,1:3,m)-r(j,1:3,k));
            % Potential between two ellipsoidal bodies - 2nd order expansion in MOI - Valid if
            % m1=m2=m3
%             Uint(k,m)=-G*mass(k)/rrel * (mass(m)+1/(2*rrel^2) * (Isum(k)+Isum(m)...
%                 -3/2*(Ivec(k,1)+Ivec(k,2)+Ivec(m,1)+Ivec(m,2)-...
%                     cos(2*phi(j,k,m))*(Ivec(k,2)-Ivec(k,1))-...
%                         cos(2*phi(j,m,k))*(Ivec(m,2)-Ivec(m,1))))); 
%             Uint(m,k)=Uint(k,m);
%             U(j)=Uint(k,m)+U(j);
            U(j)=U(j)+(-G*mass(k)*mass(m)/rrel);
        end
    end
    
    Hmag(j)=norm(H(j,1:3)); %Total system angular momentum, km^2/s
    Hrvmag(j)=norm(Hrv(j,1:3)); HwImag(j)=norm(HwI(j,1:3)); 
    E(j,1)=T(j)+U(j); %Total system energy, km^2/s^2
    
    %Energy Difference between timesteps
    if j>1
        Ediff(j-1) = E(j)-E(j-1); Hdiff(j-1) = Hmag(j)-Hmag(j-1);
    end
    
    ct=1;
    for k=1:N-1 %Compute normal directions and relative v_par & v_perp for each body
        for m=k+1:N
            rkm=r(j,:,m)-r(j,:,k);
            v_rel=v(j,:,m)-v(j,:,k);
            nhat=rkm/norm(rkm); that=cross([0 0 1], nhat);
            vpar(j,ct) = dot(v_rel,nhat);
            vperp(j,ct) = dot(v_rel,that);
            ct=ct+1;
        end
    end
end

% theta_form=zeros(length(t1),3);
% for j=1:length(t1)
%     if norm(r(j,:,1)) < 1e-10
%         theta_form(j,1) = 0; %Body 1 & 2
%         theta_form(j,2) = real(acosd(norm(dot(r(j,:,2),r(j,:,3)))/(norm(r(j,:,2))*norm(r(j,:,3))))); %Body 2 & 3
%         theta_form(j,3) = 0; %Body 3 & 1
%     else
%         theta_form(j,1) = real(acosd(norm(dot(r(j,:,1),r(j,:,2)))/(norm(r(j,:,1))*norm(r(j,:,2))))); %Body 1 & 2
%         theta_form(j,2) = real(acosd(norm(dot(r(j,:,2),r(j,:,3)))/(norm(r(j,:,2))*norm(r(j,:,3))))); %Body 2 & 3
%         theta_form(j,3) = real(acosd(norm(dot(r(j,:,3),r(j,:,1)))/(norm(r(j,:,3))*norm(r(j,:,1))))); %Body 3 & 1
%     end
% end
formCase = 'ER';
switch formCase
    case 'ER'   %Form factor is angle made between r_ab & r_ac
        a = 2; % a = center body, change based on results
        b = 1; c = 3;
        theta_form=zeros(length(t1),1);
        for j=1:length(t1)
            r_ab = r(j,:,b) - r(j,:,a);
            r_ac = r(j,:,c) - r(j,:,a);
            FF = acos(dot(r_ab,r_ac)/(norm(r_ab)*norm(r_ac))) * 180/pi;
            theta_form(j,1) = FF;
        end
    case 'LR' %Angle between Arbitrary body and other 2
        a = 1; b = 2; c = 3;
        theta_form=zeros(length(t1),1);
        for j=1:length(t1)
            r_ab = r(j,:,b) - r(j,:,a);
            r_ac = r(j,:,c) - r(j,:,a);
            FF = acos(dot(r_ab,r_ac)/(norm(r_ab)*norm(r_ac))) * 180/pi;
            theta_form(j,1) = FF;
        end
    case 'AM'
%         a = 2; %Single orbiting body
%         b = 1; c = 3;
%         theta_form=zeros(length(t1),1); FF1=zeros(length(t1),1); FF2=zeros(length(t1),1);
%         for j=1:length(t1)
%             r_ab = r(j,:,b) - r(j,:,a);
%             r_ac = r(j,:,c) - r(j,:,a);
%             FF1(j,1) = acos(dot(r(k,:,a),r_ab)/(norm(r(k,:,a))*norm(r_ab))) * 180/pi;
%             FF2(j,1) = acos(dot(r(k,:,a),r_ac)/(norm(r(k,:,a))*norm(r_ac))) * 180/pi;
%             theta_form(j,1) = FF1(j,1) + FF2(j,1);
%         end
        a = 2; b = 1; c = 3;
        theta_form=zeros(length(t1),1);
        for j=1:length(t1)
            r_ab = r(j,:,b) - r(j,:,a);
            r_ac = r(j,:,c) - r(j,:,a);
            FF = acos(dot(r_ab,r_ac)/(norm(r_ab)*norm(r_ac))) * 180/pi;
            theta_form(j,1) = FF;
        end
end

% Calculate Tidal Torque and Edot
epsilon = 5e-3; % Specific tidal energy dissipation constant
klove = 1e-5; %.0042; % Tidal love number (.03 = moon)
rho = 2e12; %[kg/km^3], density of asteroid
Gamma = zeros(length(t1),3,3); Edot = zeros(length(t1),3,3); EdotSum=zeros(length(t1),1);
for j=1:length(t1)
    for k=1:N
        for l=1:N
            if k==l
                continue
            end
            rkl = r(j,:,k)-r(j,:,l);
            Gamma(j,k,l) = sign(phi_dot(j,k,l))*3/2*klove*(3/(4*pi*rho))^2*G*mass(k)^2*mass(l)^2/...
            (norm(rkl)^6*R(k))*sin(2*epsilon); %Tidal torque
            % Must check if phi_dot is passing near 0
            delta = norm(Gamma(j,k,l)/Ivec(k,3)*sqrt(6/(pi*G*rho))); %Small Characteristic spin rate
            if abs(phi_dot(j,k,l)) < delta %Must use modified tidal torque equation
                Gamma(j,k,l)=sqrt(pi*G*rho/6)*Ivec(k,3)*phi_dot(j,k,l);
            end
            
            Edot(j,k,l) = -Gamma(j,k,l)*phi_dot(j,k,l);
            EdotSum(j,1) = EdotSum(j,1) + Edot(j,k,l);
        end
    end
end

% Calculate Hbar2, Ebar, and w_req for terminal configuration
% formCase = 'ER';
Hbar2=Hmag.^2/(G*mass(1)^3 * 2*radius);
switch formCase
    case 'ER'
        Ebar = Hbar2/(2*2.3) - 5/2;
        w_req = Hmag(end)/(3*Ivec(1,1)+2*mass(1)*(radius*2)^2);
    case 'LR'
        Ebar = Hbar2/(2*1.3) - 3;
        w_req = Hmag(end)/(3*(Ivec(1,1)+mass(1)*norm(r(end,:,1))^2));
    case 'AM'
        %Rbar must be solved for numerically based on Hbar2
        Rbar = fzero(@(Rbar) 2*(1.2+Rbar^2)^2/(3*Rbar)*(1/(Rbar-1/2)^2+1/(Rbar+1/2)^2) - Hbar2(end),3);
        Ebar = 3*Hbar2/(4*(1.2+Rbar^2)) - 1 - 1/(Rbar-1/2) - 1/(Rbar+1/2);
        [w_req, ~] = wvHtidal(Hmag(end), mass(1), Ivec(1,1), r(end,:,1), r(end,:,2), r(end,:,3));
end
EbarSI = Ebar*G*mass(1)^2;

% Calculate psi_ddot of each body based on the sum of tidal torques vs. numerical approximation
for k=1:N
    slct=[1 2 3]; slct(k)=[]; a=slct(1); b=slct(2);
    psi_ddot(1,k) = -(Gamma(end,k,a)+Gamma(end,k,b))/Ivec(k,1);
    t_tidal(1,k) = (w_req-w(end,3,k))/psi_ddot(1,k);
end
t_tidal_yrs = t_tidal/(60^2*24*365);

% Calculate psi_ddot of each body based on the sum of tidal torques using an average
% for k=1:N
%     slct=[1 2 3]; slct(k)=[]; a=slct(1); b=slct(2);
%     psi_ddot_avg(1,k) = sum(-(Gamma(tstart:tfin,k,a)+Gamma(tstart:tfin,k,b))/Ivec(k,1))/(t1(tfin)-t1(tstart));
%     t_tidal_avg(1,k) = (w_req-w(end,3,k))/psi_ddot(1,k);
% end
% t_tidal_yrs_avg = t_tidal_avg/(60^2*24*365);

psi_ddot_numer=zeros(length(t1)-1,3);
for j=1:length(t1)-1
    for k=1:3
        psi_ddot_numer(j,k) = (w(j+1,3,k)-w(j,3,k))/(t1(j+1)-t1(j));
    end
end

% Calculate system variance of Rbar_measured from the ideal Rbar (only for Mixed configurations)
if strcmp(formCase,'AM')
a = 2; %isolated body
b = 1; c = 3; %Two bodies in a contact pair
rbc_CM(:,1:3) = (r(:,:,b) + r(:,:,c))/2; %Vector pointing to CM of contact pair
Rnorm_numer(:,1) = (sqrt(sum(rbc_CM(:,:).^2,2)) + sqrt(sum(r(:,:,a).^2,2)))/(2*radius);
end

clear rrel
