%This function initilizes the configuration for the 3Body simulation.
%   config = 2 letter code for 1 of 7 possible initial configurations
%   param = The governing parameter for the configuration (either Hbar^2 or Rnorm)
%   param2 = Any secondary governing parameter (VRest angle)
%   SI = 1 to use SI units, 0 to use canonical units
%   doff = spatial displacement of one particular body from equilibrium (will speed up transition)
%   angoff = angular displacement of one particular body from tidal equilibrium (will increase
%             energy decay rate)

function [Body, G] = testConfigsBatch(input)
config = input{1}; param = input{2}; param2 = input{3}; SI = input{4}; doff = input{5}; angoff = input{6};
config = config(1:2);
global TU MU DU
G=6.67384e-20; % km^3/(kg-s^2)
rho=2e12; %kg/km^3
a=2; %radius along Body 'a' axis (stretched)
f=0; %Ellipsoidal Flattening
b=a-a*f; c=b; %radius along Body 'b' axis (compressed), c = b
mass=rho*4/3*pi*a*b*c;
radius=(3*mass/(4*rho*pi))^(1/3); %[km], Asteroid radius if perfect sphere with same mass

if SI == 0
    MU = mass; %canonical mass unit
    DU = radius; %canonical distance unit
    TU = sqrt(DU^3/(G*MU));
    G = 1;
    a=a/DU; b=b/DU; c=c/DU; radius = radius/DU; mass = mass/MU;
else
    MU =1; DU=1; TU=1;
end

for i=1:3
    Body(i).semiaxis(1:3)=[a b c];
    Body(i).R=radius;
    Body(i).mass=mass;
    Body(i).Ivec=[1/5*mass*(b^2+c^2); 1/5*mass*(a^2+c^2); 1/5*mass*(a^2+b^2)]; % MOI about ellipsoid axes
    Body(i).I=2/5*Body(i).mass*Body(i).R^2; % MOI if asteroid were perfect sphere with equivalent mass (used for initial spin/velocity calculations)
end


switch upper(config)

    case 'TM' %Transverse Mixed
% ******************** Transverse Mixed Configuration ***************************
                                % Always unstable
        Rnorm=param; %Requires physical Rnorm > sqrt(3)/2
        d=2*radius;
        Hbar2=4*(1.2+Rnorm^2)^2/(3*(Rnorm^2+.25)^(3/2)); H2=Hbar2*G*mass^3*d; Hreq=sqrt(H2);
        r1=[0,-d/2,0]; r2=[0,d/2,0]; r3=[Rnorm*d, 0, 0]; %Initial Body positions (Origin between Body 1 & 2)
        [w_req, OutputBody] = wvHtidal(Hreq, mass, Body(1).Ivec(3), r1, r2, r3);
        for k=1:length(Body)
            Body(k).r=zeros(1,3); Body(k).v=zeros(1,3); Body(k).w=zeros(1,3); Body(k).psi=zeros(1,3);
            Body(k).r=OutputBody(k).rB; %initial position, km
            Body(k).v=OutputBody(k).v_vec; %initial velocity, km/s
            Body(k).w=[0,0,w_req]; %initial angular velocity, radians/s
            Body(k).psi=0; %initial angle, radians
        end
    case 'AM' %Aligned Mixed
% ******************** Aligned Mixed Configuration ******************************
                            % Globl Emin for Rnorm>~3.09 
                            % Body 3 has offset
        Rnorm=param; %Requires physical Rnorm > 3/2
        d=2*radius;
        offset2 = angoff; offset = doff;
        Hbar2=2*(1.2+Rnorm^2)^2/(3*Rnorm)*(1/(Rnorm-1/2)^2+1/(Rnorm+1/2)^2); H2=Hbar2*G*mass^3*(2*radius); Hreq=sqrt(H2);
        r1=[-d/2, 0, 0]; r2=[d/2, 0, 0]; r3=[Rnorm*d, 0, 0]; %Initial Body positions (Origin between Body 1 & 2)
        [w_req, OutputBody] = wvHtidal((1-offset)*Hreq, mass, Body(1).I, r1, r2, r3);
        for ii=1:length(Body)
            Body(ii).r(1,1:3)=OutputBody(ii).rB; %initial position, km
            Body(ii).v(1,1:3)=OutputBody(ii).v_vec; %initial velocity, km/s
            Body(ii).w(1,1:3)=[0,0,w_req]; %initial angular velocity, radians/s
            Body(ii).psi(1,1:3)=[0,0,0]; %initial angle, radians
        end
        Body(3).psi(1,1:3) = Body(3).psi(1,1:3) + offset2;
    case 'LO' %Lagrange Orbiting
%                 ********** Lagrange Orbiting (Triangle) **********
        Rnorm=param; Hbar2=(3*Rnorm^2+.9)^2/(3*Rnorm^3); H2=Hbar2*G*mass^3*(2*radius); Hreq=sqrt(H2);
        d=2*radius*Rnorm;
        Body(1).r(1,1:3)=d*[ -1/2, -sqrt(3)/6, 0]; %initial position, km
        Body(2).r(1,1:3)=d*[ 0, sqrt(3)/3, 0];
        Body(3).r(1,1:3)=d*[ 1/2, -sqrt(3)/6, 0];

        I1cm=Body(1).Ivec(3)+mass*norm(Body(1).r)^2; %Moment of Ineteria of Body about CM (Parallel axis Thm)
        w_req=Hreq/(3*I1cm); v_req=w_req*norm(Body(1).r);
        %vt=sqrt(norm(Body(1).r)*(2*G*mass*cosd(30)/d^2)); wsys=vt/norm(Body(1).r); %Force balance, Alt. Method
        wvec=w_req*[0,0,1];

        Body(1).v(1,1:3)=v_req*cross(wvec,Body(1).r)/norm(cross(wvec,Body(1).r));
        Body(2).v(1,1:3)=v_req*cross(wvec,Body(2).r)/norm(cross(wvec,Body(2).r));
        Body(3).v(1,1:3)=v_req*cross(wvec,Body(3).r)/norm(cross(wvec,Body(3).r));
        Body(1).w(1,1:3)=[0,0,w_req]; %initial angular velocity, radians/s
        Body(1).psi=pi/6; %initial bulge angle, radians
        Body(2).w(1,1:3)=[0,0,w_req]; %initial angular velocity, radians/s
        Body(2).psi=-pi/2; %initial bulge angle, radians
        Body(3).w(1,1:3)=[0,0,w_req]; %initial angular velocity, radians/s
        Body(3).psi=5*pi/6; %initial bulge angle, radians
    case 'EO' %Euler Orbiting
%               ********** Euler Orbiting (Line) - Always unstable ********** 
        Rnorm=param; %normalized separation distance between inner and outer bodies, R>=1 
        Hbar2=5/36*(6*Rnorm^2+.9)^2/Rnorm^3; H2=Hbar2*G*mass^3*(2*radius); Hreq=sqrt(H2);
        d=2*radius*Rnorm;

        w_req=Hreq/(3*Body(1).Ivec(3)+2*mass*d^2); %Required spin rate for system to be stable (using Parallel axis Thm)
        v_req=w_req*d; %Required velocity mag of adjacent bodies to spin at correct rate
        % vt=sqrt(d*G*mass*(1/d^2+1/(2*d)^2)); wsys=vt/d; % Using force balance, to validate as Alt. Method
        
        Body(1).r=[0,0,0]; Body(2).r=[d,0,0]; Body(3).r=[-d,0,0];
        Body(1).v=[0,0,0]; Body(2).v=[0,v_req,0]; Body(3).v=[0,-v_req,0];
        Body(1).w=[0,0,w_req]; Body(2).w=[0,0,w_req]; Body(3).w=[0,0,w_req];
        Body(1).psi=0; Body(2).psi=0; Body(3).psi=0;
    case 'LR' %Lagrange Resting
%                 ********** Lagrange Resting (Triangle) **********
                        % Global Emin for H^2 < 1.98375
                        % Body 3 has offset
        d=2*radius;
        offset = d*doff; offset2 = angoff;
        Hbar2= param; 
        H2=Hbar2*G*mass^3*d; Hreq=sqrt(H2);
        Body(1).r(1,1:3)=d*[ -1/2, -sqrt(3)/6, 0]; %initial position, km
        Body(2).r(1,1:3)=d*[ 0, sqrt(3)/3, 0];
        Body(3).r(1,1:3)=d*[ 1/2, -sqrt(3)/6, 0];

        I1cm=Body(1).Ivec(3)+mass*norm(Body(1).r)^2; %Moment of Inertia about CM (Parallel axis Thm)
        w_req=Hreq/(3*I1cm); v_req=w_req*norm(Body(1).r);
        wvec=w_req*[0 0 1]';
        if Hbar2==0
            Body(1).v(1,1:3)=[0,0,0]; Body(2).v(1,1:3)=[0,0,0]; Body(3).v(1,1:3)=[0,0,0];
        else
            Body(1).v(1,1:3)=v_req*cross(wvec,Body(1).r)/norm(cross(wvec,Body(1).r));
            Body(2).v(1,1:3)=v_req*cross(wvec,Body(2).r)/norm(cross(wvec,Body(2).r));
            Body(3).v(1,1:3)=v_req*cross(wvec,Body(3).r)/norm(cross(wvec,Body(3).r));
        end
        Body(1).w(1,1:3)=[0,0,w_req]; %initial angular velocity, radians/s
        Body(1).psi(1,1:3)=pi/6; %initial angle, radians
        Body(2).w(1,1:3)=[0,0,w_req]; %initial angular velocity, radians/s
        Body(2).psi(1,1:3)=-pi/2; %initial angle, radians
        Body(3).w(1,1:3)=[0,0,w_req]; %initial angular velocity, radians/s
        Body(3).psi(1,1:3)=5*pi/6 + offset2; %initial angle, radians
        Body(3).r(1,1)=Body(3).r(1,1) + offset;
    case 'ER' %Euler Resting
%                 ********** Euler Resting (Line) **********
                    % Global Emin for 2.99< H^2 < 5.65907
                    % Body 2 has offset
        d=2*radius;
        offset = d*doff; offset2 = angoff;
        Hbar2= param;
        H2=Hbar2*G*mass^3*(d); Hreq=sqrt(H2);
        w_req=Hreq/(3*Body(1).Ivec(3)+2*mass*d^2); %Required spin rate for system to be stable
        v_req=w_req*d; %Required velocity mag of adjacent bodies to spin at correct rate
        Body(1).r=[0,0,0]; Body(2).r=[d+offset,0,0]; Body(3).r=[-d,0,0];
        Body(1).v=[0,0,0]; Body(2).v=[0,v_req,0]; Body(3).v=[0,-v_req,0];
        Body(1).w=[0,0,w_req]; Body(2).w=[0,0,w_req]; Body(3).w=[0,0,w_req];
        Body(1).psi=0; Body(2).psi=offset2; Body(3).psi=0;
    case 'VR' %V Resting
%                 ********** V Resting (V Shape) **********
                            % Unstable
        restAngle = param2; %Angle 3rd body makes with other 2 bodies
        Hbar2 = param;
%         Hbar2=(2.9+4*sind(restAngle/2)^2)^2/(24*sind(restAngle/2)^3);
        H2=Hbar2*G*mass^3*(2*radius); Hreq=sqrt(H2);
        d=2*radius;

        Body(1).r=[-d,0,0]; Body(2).r=[0,0,0]; Body(3).r=[-d*cosd(restAngle), d*sind(restAngle), 0]; %Initial Body positions (Origin at body 2)
        RCM=(Body(1).r+Body(2).r+Body(3).r)/3; % Barycenter for system
        Body(1).r=Body(1).r-RCM; Body(2).r=Body(2).r-RCM; Body(3).r=Body(3).r-RCM; %Convert to Barycenter coordinates

        w_req=Hreq/(3*Body(1).Ivec(3)+mass*(norm(Body(1).r)^2+norm(Body(2).r)^2+norm(Body(3).r)^2)); %Required spin rate of bodies
        v1=w_req*norm(Body(1).r); v2=w_req*norm(Body(2).r); v3=w_req*norm(Body(3).r);  %Required tangential velocities
        w_vec=[0,0,1]; %Rotation axis
        vvec_1=cross(w_vec,Body(1).r)/norm(cross(w_vec,Body(1).r)); %Direction of veloicty vector for body 1
        vvec_2=cross(w_vec,Body(2).r)/norm(cross(w_vec,Body(2).r)); %Direction of veloicty vector for body 2
        vvec_3=cross(w_vec,Body(3).r)/norm(cross(w_vec,Body(3).r)); %Direction of veloicty vector for body 3
        Body(1).v=v1*vvec_1; Body(2).v=v2*vvec_2; Body(3).v=v3*vvec_3;
        Body(1).w=[0,0,w_req]; Body(2).w=[0,0,w_req]; Body(3).w=[0,0,w_req];
        Body(1).psi=[0,0,0]; Body(2).psi=[0,0,0]; Body(3).psi=[0,0,0];

end


