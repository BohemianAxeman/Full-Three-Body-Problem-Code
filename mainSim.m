
%%%%%%%%%%%% Orbit Propagation Using Custom RK45 Integrator %%%%%%%%%%%%
global MU DU TU
tend=1*31536000/TU; %31536000;
step= 79/TU; %initial step size to try (or to use for fixed step propagation)
ColorSpace=distinguishable_colors(N+1);
global currTime
tol=1e-6;
varstep=0; %0=fixed step, 1=variable step

%Initialize contact and collision property
for a=1:N
    Body(a).contact(1,1:N)=0;
    Body(a).collision(1,1:N)=0;
end

tinit = clock;
% [tvec, BodyOut] = mainPropagator(tend,Body,G,varstep,step,tol); % RK variable timestep
[tvec, BodyOut] = mainPropagatorVerlet(tend,Body,G,varstep,step,tol); % Verlet
% [tvec, BodyOut] = SHAKE_prop(tend,Body,G,varstep,step,tol)
tfinal = clock;
telap = etime(tfinal,tinit);
disp(['Elapsed time: ', num2str(telap),' seconds'])

processData
% plotData
