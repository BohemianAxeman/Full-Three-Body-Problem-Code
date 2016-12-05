
%%%%%%%%%%%%%%%% Data Visualization %%%%%%%%%%%%%%%% 
% Must have run processData first
co = [    0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];
set(groot,'defaultAxesColorOrder',co)

global TU MU DU
G=6.67384e-20; % km^3/(kg-s^2)
N = 3;
ColorSpace=co;


%Plot contact condition
figure
for k=1:N-1
    Ba=mod(k,N)+1; Bb=mod(k+1,N)+1;
    subplot(N-1,N-1,2*k-1)
    plot(t1,contacts(:,Ba,k)), title(['B',num2str(k),' & B',num2str(Ba)])
    subplot(N-1,N-1,2*k)
    plot(t1,contacts(:,Bb,k)), title(['B',num2str(k),' & B',num2str(Bb)])
end

%Plot orbit positions
figure
hold all
for k=1:N
    plot3(r(:,1,k),r(:,2,k),r(:,3,k),'.-','LineWidth',1,'MarkerSize',5,'Color',ColorSpace(k,:))
end
%plot3(Rcm(:,1),Rcm(:,2),Rcm(:,3),'kx')
xlabel('x position')
ylabel('y position')
zlabel('z position')
legend('Primary Body','Secondary Body','Tertiary Body')
title('Generic 3 Body System')
axis equal

    %Plot energy and angular momentum of system
figure
subplot(2,2,1)
plot(t1,E)
xlabel('Time (sec)')
ylabel('Energy (kg*km^2/s^2)')
title('System Energy')
%xlim([t1(1) t1(end)]); ylim([min(E) max(E)]);
grid on
subplot(2,2,2)
plot(t1,E-E(1))
xlabel('Time (sec)')
ylabel('Energy (kg*km^2/s^2)')
title('Relative System Energy (E-E0)')
%xlim([t1(1) t1(end)]); ylim([min(E-E(1))-.01*abs(min(E-E(1))) max(E-E(1))+.01*abs(max(E-E(1)))]);
grid on
subplot(2,2,3)
plot(t1,Hmag)
xlabel('Time (sec)')
ylabel('Angular Momentum (kg*km^2/s)')
title('System Angular Momentum')
xlim([t1(1) t1(end)]); %ylim([min(Hmag)-.01*abs(min(Hmag)) max(Hmag)+.01*abs(max(Hmag))]);
grid on
subplot(2,2,4)
plot(t1,Hmag-Hmag(1))
xlabel('Time (sec)')
ylabel('Angular Momentum (kg*km^2/s)')
title('Relative System Angular Momentum (H-H0)')
xlim([t1(1) t1(end)]);  %ylim([min(Hmag-Hmag(1))-.01*abs(min(Hmag-Hmag(1))) max(Hmag-Hmag(1))+.01*abs(max(Hmag-Hmag(1)))]);
grid on

        % Plot Just Relative Energy & AM
        tfin=find(t1>=1.2e5,1);
figure
subplot(1,2,1)
plot(t1,E)
xlabel('Time (sec)')
ylabel('Energy (kg*km^2/s^2)')
title('System Energy')
grid on
xlim([0 t1(tfin)])
subplot(1,2,2)
plot(t1,E-E(1))
xlabel('Time (sec)')
ylabel('Energy (kg*km^2/s^2)')
title('Relative System Energy (E-E0)')
grid on
xlim([0 t1(tfin)])
% subplot(1,2,2)
% plot(t1,Hmag-Hmag(1))
% xlabel('Time (sec)')
% ylabel('Energy (kg*km^2/s^)')
% title('Relative System AM (H-H0)')


%     %Plot energy and angular momentum of system for a specific range of times
% tstart=find(t1>=1/2*t1(end),1); tfin=find(t1>=t1(end),1);
tstart=find(t1>=0,1); tfin=find(t1>=t1(end),1);
figure
subplot(2,1,1)
plot(t1(tstart:tfin),(E(tstart:tfin)-E(tstart)))
xlabel('Time (sec)')
ylabel('Energy (kg*km^2/s^2)')
title('Error in System Energy'); grid on
% axis([t1(tstart),t1(tfin),-.1,.15])
subplot(2,1,2)
plot(t1(tstart:tfin),Hmag(tstart:tfin)-Hmag(tstart))
xlabel('Time (sec)')
ylabel('Angular Momentum (kg*km^2/s)')
title('Error in System Angular Momentum'); grid on
% axis([t1(tstart),t1(tfin),-1000,2000])

    %Plot Absolute and Relative Kinetic & Potential Energy
%Absolute Energy
    figure
    subplot(3,2,1)
    plot(t1,T)
    xlabel('Time (sec)')
    ylabel('T (kg*km^2/s^2)')
    title('Kinetic Energy'), grid on
    subplot(3,2,3)
    plot(t1,U)
    xlabel('Time (sec)')
    ylabel('U (kg*km^2/s^2)')
    title('Potential Energy'), grid on
    subplot(3,2,5)
    plot(t1,E)
    xlabel('Time (sec)')
    ylabel('E (kg*km^2/s^2)')
    title('Total System Energy'), grid on
%Relative Energy
    subplot(3,2,2)
    plot(t1,T-T(1))
    xlabel('Time (sec)')
    ylabel('K-K0 (kg*km^2/s^2)')
    title('Relative Kinetic Energy'), grid on
    subplot(3,2,4)
    plot(t1,U-U(1))
    xlabel('Time (sec)')
    ylabel('U-U0 (kg*km^2/s^2)')
    title('Relative Potential Energy'), grid on
    subplot(3,2,6)
    plot(t1,E-E(1))
    xlabel('Time (sec)')
    ylabel('E-E0 (kg*km^2/s^2)')
    title('Relative Total System Energy'), grid on
    
    
    %Plot Components of Kinetic Energy & Angular Momentum
%Kinetic Energy
    figure
    subplot(3,2,1)
    plot(t1,Ttrans)
    xlabel('Time (sec)')
    ylabel('T_{trans} (kg*km^2/s^2)')
    title('Translational Kinetic Energy')
    subplot(3,2,3)
    plot(t1,Trot)
    xlabel('Time (sec)')
    ylabel('T_{rot} (kg*km^2/s^2)')
    title('Rotational Kinetic Energy')
    subplot(3,2,5)
    plot(t1,T)
    xlabel('Time (sec)')
    ylabel('T (kg*km^2/s^2)')
    title('Total Kinetic Energy')
%Angular Momentum
    subplot(3,2,2)
    plot(t1,Hrvmag)
    xlabel('Time (sec)')
    ylabel('AM (kg*km^2/s)')
    title('Rotation About Barycenter (r x v)')
    subplot(3,2,4)
    plot(t1,HwImag)
    xlabel('Time (sec)')
    ylabel('AM (kg*km^2/s)')
    title('Rotation About Body Axis (\psi * I)')
    subplot(3,2,6)
    plot(t1,Hmag)
    xlabel('Time (sec)')
    ylabel('AM (kg*km^2/s)')
    title('Total Angular Momentum')
    
    %Plot Components of Kinetic Energy & Angular Momentum
%Kinetic Energy
    figure
    subplot(3,2,1)
    plot(t1,Ttrans)
    xlabel('Time (sec)')
    ylabel('T_{trans} (kg*km^2/s^2)')
    title('Translational Kinetic Energy')
    subplot(3,2,3)
    plot(t1,Trot)
    xlabel('Time (sec)')
    ylabel('T_{rot} (kg*km^2/s^2)')
    title('Rotational Kinetic Energy')
    subplot(3,2,5)
    plot(t1,T)
    xlabel('Time (sec)')
    ylabel('T (kg*km^2/s^2)')
    title('Total Kinetic Energy')
%Angular Momentum
    subplot(3,2,2)
    plot(t1,Hrvmag)
    xlabel('Time (sec)')
    ylabel('AM (kg*km^2/s)')
    title('Rotation About Barycenter (r x v)')
    subplot(3,2,4)
    plot(t1,HwImag)
    xlabel('Time (sec)')
    ylabel('AM (kg*km^2/s)')
    title('Rotation About Body Axis (\psi * I)')
    subplot(3,2,6)
    plot(t1,Hmag)
    xlabel('Time (sec)')
    ylabel('AM (kg*km^2/s)')
    title('Total Angular Momentum')
    
     %Plot Relative Components of Angular Momentum for specific range of time
%Angular Momentum
%     tstart=find(t1>=1/2*t1(end),1); tfin=find(t1>=t1(end),1);
    figure
    subplot(3,1,1)
    plot(t1(tstart:tfin),Hrvmag(tstart:tfin)-Hrvmag(tstart))
%     xlabel('Time (sec)')
    ylabel('AM (kg*km^2/s)')
    title('Relative AM (r x v)'); grid on
    subplot(3,1,2)
    plot(t1(tstart:tfin),HwImag(tstart:tfin)-HwImag(tstart))
%     xlabel('Time (sec)')
    ylabel('AM (kg*km^2/s)')
    title('Relative AM (\psi * I)'); grid on
    subplot(3,1,3)
    plot(t1(tstart:tfin),Hmag(tstart:tfin)-Hmag(tstart))
    xlabel('Time (sec)')
    ylabel('AM (kg*km^2/s)')
    title('Relative Total AM'); grid on
    
%     %Plot AM (r x v) for each body
%     figure
%     for i=1:N
%        subplot(1,3,i)
%        plot(t1,Hrv_multi(:,3,i))
%        xlabel('Time (sec)')
%        ylabel('AM (kg*km^2/s)')
%        title(['Body ',num2str(i),' AM (r x v)'])
%     end
    
    %Plot CM-Relative Components of Absolute & Relative Kinetic Energy & Angular Momentum
%Kinetic Energy
    figure
    subplot(2,2,1)
    plot(t1,Trel_CM)
    xlabel('Time (sec)')
    ylabel('T_{trans} (kg*km^2/s^2)')
    title('Translational Kinetic Energy (v-vCM)')
    subplot(2,2,3)
    plot(t1,Trel_CM-Trel_CM(1))
    xlabel('Time (sec)')
    ylabel('T_{trans} (kg*km^2/s^2)')
    title('Relative Translational Kinetic Energy (v-CM)-T0')
%Angular Momentum
    subplot(2,2,2)
    plot(t1,sqrt(sum(Hrel_CM.^2,2)))
    xlabel('Time (sec)')
    ylabel('AM (kg*km^2/s)')
    title('Rotation About Barycenter (r-Rcm x v-vCM)')
    subplot(2,2,4)
    plot(t1,sqrt(sum(Hrel_CM.^2,2))-sqrt(sum(Hrel_CM(1,:).^2,2)))
    xlabel('Time (sec)')
    ylabel('AM (kg*km^2/s)')
    title('Relative Rotation About Barycenter (r-Rcm x v-vCM)-H0')
    
    % Plot energy gain between Time steps
    figure
    Ediff_abs=abs(Ediff); Hdiff_abs=abs(Hdiff);
    subplot(2,2,1)
    semilogy(t1(2:end),Ediff_abs); grid on
    xlabel('Time (sec)'); ylabel('\Delta E (kg*km^2/s^2)'); title('Energy Diff. versus Time')
    subplot(2,2,2)
    semilogy([2:length(t1)],Ediff_abs); grid on
    xlabel('iteration #'); ylabel('\Delta E (kg*km^2/s^2)'); title('Energy Diff. versus Iteration')
    subplot(2,2,3)
    semilogy(t1(2:end),Hdiff_abs); grid on
    xlabel('Time (sec)'); ylabel('\Delta AM (kg*km^2/s)'); title('AM Diff. versus Time')
    subplot(2,2,4)
    semilogy([2:length(t1)],Hdiff_abs); grid on
    xlabel('iteration #'); ylabel('\Delta AM (kg*km^2/s)'); title('AM Diff. versus Iteration')
    
    % Plot energy gain between Time steps
    figure
    Ediff_abs=abs(Ediff); Hdiff_abs=abs(Hdiff);
    subplot(1,2,1)
    semilogy(t1(2:end),Ediff_abs); grid on
    xlim([0 t1(end-1)])
    xlabel('Time (sec)'); ylabel('\Delta E (kg*km^2/s^2)'); title('Energy Diff. versus Time')
    subplot(1,2,2)
    semilogy(t1(2:end),Hdiff_abs); grid on
    xlim([0 t1(end-1)])
    xlabel('Time (sec)'); ylabel('\Delta AM (kg*km^2/s)'); title('AM Diff. versus Time')

    %Plot position of Body axis and Body axis spin rate
figure
for i=1:N
    subplot(2,N,i)
    plot(t1,psi(:,1,i))
    hold all
    xlabel('Time (sec)')
    ylabel('Body Axis Angle \psi [rad]')
    title(['Body: ',num2str(i)])
    grid on
    subplot(2,N,i+N)
    plot(t1,w(:,3,i))
    xlabel('Time (sec)')
    ylabel('Body Axis Spin Rate $\dot{\psi}$ [rad/s]','interpreter','latex')
    grid on
    %axis([t1(1) t1(end) -2*max(w(:,3,i)) 2*max(w(:,3,i))])
end

    %Plot Relative Body Spin Angle and Spin Rate, Phi & Phi_dot
figure
for i=1:3
    slct=1:3; slct(i)=[];
    
    subplot(2,6,2*i-1)
    plot(t1,phi_dot(:,i,slct(1)))
    hold all
    xlabel('Time (s)')
    if i==1
        ylabel('Relative Spin Rate $\dot{\phi}$ [rad/s]','interpreter','latex')
    end
    xlim([t1(1) t1(end-1)]); %ylim([min(phi_dot(:,i,slct(1))) max(phi_dot(:,i,slct(1)))]);
    title(['$\dot{\phi}_{',num2str(i),num2str(slct(1)),'}$'],'interpreter','latex')
    grid on
    subplot(2,6,2*i)
    plot(t1,phi_dot(:,i,slct(2)))
    xlabel('Time (s)')
    xlim([t1(1) t1(end-1)]); %ylim([min(phi_dot(:,i,slct(2))) max(phi_dot(:,i,slct(2)))]);
    title(['$\dot{\phi}_{',num2str(i),num2str(slct(2)),'}$'],'interpreter','latex')
    grid on
    
    subplot(2,6,2*i+5)
    plot(t1,phi(:,i,slct(1)))
    hold all
    xlabel('Time (s)')
    ylabel('Relative Body Angle \phi [rad]')
    grid on
    subplot(2,6,2*i+6)
    plot(t1,phi(:,i,slct(2)))
    xlabel('Time (s)')
    ylabel('Relative Body Angle \phi [rad]')
    grid on
end

 %Plot Relative Body Spin Angle and Spin Rate, Phi & Phi_dot
figure
subplot(1,3,1)
plot(t1,phi_dot(:,1,2))
xlabel('Time (s)')
ylabel('Relative Spin Rate $\dot{\phi}$ [rad/s]','interpreter','latex')
xlim([t1(1) t1(end-1)]); %ylim([min(phi_dot(:,i,slct(1))) max(phi_dot(:,i,slct(1)))]);
title(['$\dot{\phi}_{',num2str(1),num2str(2),'}$'],'interpreter','latex')
grid on
subplot(1,3,2)
plot(t1,phi_dot(:,2,3))
xlabel('Time (s)')
ylabel('Relative Spin Rate $\dot{\phi}$ [rad/s]','interpreter','latex')
xlim([t1(1) t1(end-1)]); %ylim([min(phi_dot(:,i,slct(1))) max(phi_dot(:,i,slct(1)))]);
title(['$\dot{\phi}_{',num2str(2),num2str(3),'}$'],'interpreter','latex')
grid on
subplot(1,3,3)
plot(t1,phi_dot(:,3,1))
xlabel('Time (s)')
ylabel('Relative Spin Rate $\dot{\phi}$ [rad/s]','interpreter','latex')
xlim([t1(1) t1(end-1)]); %ylim([min(phi_dot(:,i,slct(1))) max(phi_dot(:,i,slct(1)))]);
title(['$\dot{\phi}_{',num2str(3),num2str(1),'}$'],'interpreter','latex')
grid on




    %Plot Relative Body Axis Spin Angle and Spin Rate, Theta & Theta_dot
figure
for i=1:3
    slct=1:3; slct(i)=[];
    subplot(2,6,2*i-1)
    plot(t1,theta(:,i,slct(1)))
    hold all
    xlabel('Time (sec)')
    ylabel('Relative Axis Angle \theta [rad]')
    title(['Body: ',num2str(i),' & ',num2str(slct(1))])
    grid on
    subplot(2,6,2*i)
    plot(t1,theta(:,i,slct(2)))
    xlabel('Time (sec)')
    ylabel('Relative Axis Angle \theta [rad]')
    title(['Body: ',num2str(i),' & ',num2str(slct(2))])
    grid on
    subplot(2,6,2*i+5)
    plot(t1,theta_dot(:,i,slct(1)))
    hold all
    xlabel('Time (sec)')
    ylabel('Spin Rate $\dot{\theta}$ [rad/s]','interpreter','latex')
    grid on
    subplot(2,6,2*i+6)
    plot(t1,theta_dot(:,i,slct(2)))
    xlabel('Time (sec)')
    ylabel('Spin Rate $\dot{\theta}$ [rad/s]','interpreter','latex')
    grid on
end

%     %Plot Bodies in rotating B1 Frame
% figure
% for k=1:N
%     plot(rB1(:,1,k),rB1(:,2,k),'.-','LineWidth',1,'MarkerSize',5,'Color',ColorSpace(k,:))
%     hold all
% end
% xlabel('x position')
% ylabel('y position')
% zlabel('z position')
% legend('Primary Body','Secondary Body','Tertiary Body')
% title('Body motion in B1 frame')
% axis equal

    %Plot Distance from Bodies to Center of Mass
dCM=zeros(length(t1),N); %distance from each body to CM
figure
for k=1:N
    dCM(:,k)=sqrt(sum(abs(r(:,1:3,k)-Rcm(:,1:3)).^2,2));
    subplot(N,1,k)
    plot(t1,dCM(:,k),'Color',ColorSpace(k,:))
    xlabel('Time (sec)')
    ylabel('Distance (km)')
    title(strcat('Distance from CM to Body - ', num2str(k)))
end

%Plot of Center of Mass
figure
subplot(2,1,1)
plot(t1,RcmMAG-RcmMAG(1))
xlabel('Time (sec)')
ylabel('Distance (km)')
title('Translation of Center of Mass')
subplot(2,1,2)
plot(Rcm(:,1)-Rcm(1,1),Rcm(:,2)-Rcm(1,2))
xlabel('X (km)')
ylabel('Y (km)')
title('XY Translation of CM')

%Plot Velocities of each Body
figure
for k=1:N
    subplot(N,1,k)
    plot(t1,vMAG(:,k),'.-','Color',ColorSpace(k,:))
    xlabel('Time (sec)')
    ylabel('Velocity (km/s)')
    title(strcat('Velocity Magnitude of Body ', num2str(k)))
end

%Plot Relative Velocities of each Body
figure
z=0;
for k=1:N-1
    for j=k+1:N
        z=z+1; 
        vrel(:,z)=sqrt(sum((v(:,1:3,k)-v(:,1:3,j)).^2,2));
        subplot(N,1,z)
        plot(t1,vrel(:,z),'.-','Color',ColorSpace(z,:))
        xlabel('Time (sec)')
        ylabel('Relative Velocity (km/s)')
        title(strcat('RelVel of Bodies: ', num2str(k), '& ', num2str(j)))
    end
end
    
%Plot Relative Positions of each Body from initial starting point
figure
z=0;
for k=1:N-1
    for j=k+1:N
        z=z+1; 
        rrel(:,z)=sqrt(sum((r(:,1:3,k)-r(:,1:3,j)).^2,2));
        subplot(N,1,z)
        plot(t1,rrel(:,z)-rrel(1,z),'.-','Color',ColorSpace(z,:))
        xlabel('Time (sec)')
        ylabel('Mutual Distance (km)')
        title(['Rel. Pos. of Bodies: ', num2str(k), ' & ', num2str(j),  ' from t_0'])
    end
end

%Plot Displacement of each Body from initial starting point
figure
z=0;
for k=1:N
    z=z+1; 
    rrel(:,z)=sqrt(sum((r(:,1:3,k)).^2,2));
    subplot(N,1,z)
    plot(t1,rrel(:,z)-rrel(1,z),'.-','Color',ColorSpace(z,:))
    xlabel('Time (sec)')
    ylabel('Displacement (km)')
    title(['Displacement of Body ', num2str(k),' from t0'])
end

%Plot Relative Positions of each Body
figure
z=0;
for k=1:N-1
    for j=k+1:N
        z=z+1; 
        rrel(:,z)=sqrt(sum((r(:,1:3,k)-r(:,1:3,j)).^2,2));
        subplot(N,1,z)
        plot(t1,rrel(:,z),'.-','Color',ColorSpace(z,:))
        xlabel('Time (sec)')
        ylabel('Mutual Distance (km)')
        title(strcat('RelPos of Bodies: ', num2str(k), '& ', num2str(j)))
        grid on
    end
end



%   % Plot Relative Positions of each Body for a specific range of time
tstart=find(t1>=1/2*t1(end),1); tfin=find(t1>=t1(end),1);
figure
z=0;
for k=1:N-1
    for j=k+1:N
        z=z+1; 
        rrel(:,z)=sqrt(sum((r(:,1:3,k)-r(:,1:3,j)).^2,2));
        subplot(N,1,z)
        plot(t1(tstart:tfin),rrel(tstart:tfin,z)-rrel(tstart,z),'.-','Color',ColorSpace(z,:))
        xlabel('Time (sec)')
        ylabel('Mutual Distance (km)')
        title(strcat('RelPos of Bodies: ', num2str(k), '& ', num2str(j)))
        grid on
    end
end

%Velocity Profile of Body 1
figure
plot(v(:,1,1),v(:,2,1))
hold all
plot(v(1,1,1),v(1,2,1),'g*')
plot(v(end,1,1),v(end,2,1),'r*')
xlabel('vx (km/s)')
ylabel('vy (km/s)')
title('Velocity profile of body 1')

%Plot Time steps (for variable Time steps)
for k=1:length(t1)-1
   tstep(k,1)=t1(k+1)-t1(k); 
end
figure
subplot(2,1,1)
plot(1:length(tstep), tstep)
xlabel('iteration #')
ylabel('Time step (s)')
title('Time steps taken')
subplot(2,1,2)
plot(t1(1:end-1),tstep)
xlabel('Time (s)')
ylabel('Time step (s)')
title('Time step taken at current Time')

% Plot normal/perpendicular components of relative body velocities
figure
for k=1:N
    switch k
        case 1; text='1 & 2'; 
        case 2; text='1 & 3';     
        case 3; text='2 & 3'; 
    end
    subplot(3,2,2*k-1)
    plot(t1,vpar(:,k)) %Parallel to normal
    xlabel('Time (sec)')
    ylabel('velocity (km/s)')
    title(['Relative normal velocity, ',text])
    subplot(3,2,2*k)
    plot(t1,vperp(:,k)) %Perpendicular to normal
    xlabel('Time (sec)')
    ylabel('velocity (km/s)')
    title(['Relative tangential velocity, ',text])
end

% Plot orbital rate of bodies about barycenter
figure
for k=1:N
   nCM(:,k) = ((r(:,1,k)-Rcm(:,1)).*v(:,2,k)-(r(:,2,k)-Rcm(:,2)).*v(:,1,k))./((r(:,1,k)-Rcm(:,1)).^2 + (r(:,2,k)-Rcm(:,2)).^2);
   subplot(1,3,k)
   plot(t1,nCM(:,k))
   xlabel('Time (sec)'); ylabel('n (rad/s)'); title(['Body ',num2str(k),' mean motion about CM'])
end

% Plot Configuration Form Factor against ideal configuration
formCase = 'ER';
figure
switch formCase
    case 'LR'
        form_case = 60;
        plot(t1,theta_form(:,1))
        xlabel('Time (sec)'); ylabel('Form Factor (deg)'); grid on
        title('Form Factor Angle over time')
        legend('F\theta')
    case 'ER'
        form_case = 0;
        plot(t1,theta_form(:,1))
        xlabel('Time (sec)'); ylabel('Form Factor (deg)'); grid on
        title('Form Factor Angle over time')
        legend('F\theta')
    case 'AM'
        form_case = 0;
%         plot(t1,FF1(:,1)); hold all; 
%         plot(t1,FF2(:,1))
%         xlabel('Time (sec)'); ylabel('Form Factor (deg)'); grid on
%         title('Form Factor Angle over time')
%         legend('F\theta_{1}', 'F\theta_{3}')
        plot(t1,theta_form(:,1));
        xlabel('Time (sec)'); ylabel('Form Factor (deg)'); grid on
        title('Form Factor Angle over time')
        legend('F\theta')
end
% plot(t1,theta_form(:,1)-form_case(1), t1, theta_form(:,2)-form_case(2), t1, theta_form(:,3)-form_case(3))
% xlabel('Time (sec)'); ylabel('Form Factor (deg)'); grid on
% title('Form Factor Angle over time')
% legend('F\theta_{12}','F\theta_{23}','F\theta_{32}')

    % Plot Tidal Torque and Net Energy Loss from each Tidal Torque
figure
for i=1:3
    slct=1:3; slct(i)=[];
    
    subplot(2,6,2*i-1)
    plot(t1,Gamma(:,i,slct(1)))
    hold all
    xlabel('Time (s)')
    if i==1
        ylabel('Tidal Torque $\Gamma$','interpreter','latex')
    end
    xlim([t1(1) t1(end)]); %ylim([min(phi_dot(:,i,slct(1))) max(phi_dot(:,i,slct(1)))]);
    title(['$\Gamma_{',num2str(i),num2str(slct(1)),'}$'],'interpreter','latex')
    grid on
    subplot(2,6,2*i)
    plot(t1,Gamma(:,i,slct(2)))
    xlabel('Time (s)')
    xlim([t1(1) t1(end)]); %ylim([min(phi_dot(:,i,slct(2))) max(phi_dot(:,i,slct(2)))]);
    title(['Tidal Torque $\Gamma_{',num2str(i),num2str(slct(2)),'}$'],'interpreter','latex')
    grid on
    
    subplot(2,6,2*i+5)
    plot(t1,Edot(:,i,slct(1))*1000)
    hold all
    xlabel('Time (s)')
    ylabel('Energy Decay Rate (kW)')
    grid on
    subplot(2,6,2*i+6)
    plot(t1,Edot(:,i,slct(2))*1000)
    xlabel('Time (s)')
    ylabel('Energy Decay Rate (kW)')
    grid on
end

    % Plot Total Energy Decay Rate over time
figure
plot(t1,EdotSum*10^6)
xlabel('Time (sec)')
ylabel('Energy Decay Rate (W)')
title('Total Analytical Rate of Energy loss due to Tidal Torques')
grid on
    
%     % Plot System Energy with Minimum Energy
% figure
% plot(t1(tstart:tfin),E(tstart:tfin)); hold all
% plot([t1(tstart) t1(tend)],[EbarSI(end) EbarSI(end)])

    % Plot numerical psi_ddot
figure; hold all
for k=1:3
    plot(t1(1:end-1),psi_ddot_numer(:,k))
end

    % Plot phase space of bodies 1 & 2
figure
for k=2:3
    subplot(2,1,k-1)
    plot(r(:,1,k),v(:,1,k))
    title('Phase Space of')
end


