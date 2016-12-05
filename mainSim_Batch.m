
%%%%%%%%%%%% Orbit Propagation Using Custom Integrator %%%%%%%%%%%%
global MU DU TU currTime

% //////////    List of simulations to run in batches   //////////
%TestBatches(n) =     {1,       2,      3,      4,  5,      6,      7,      8,          9,      10};
%TestBatches(n) =     {config,  param,  param2, SI, doff,   angoff, Dura,   initStep,   tol,    varstep};
    
    % Aligned Mixed --> Euler Resting
% TestBatches(1).set =  {'AM2ER', 2.33696,    [], 0,  1e-2,   0,     2e6,    1e-3,       1e-6,   1};
% TestBatches(1).set =  {'AM2ER', 3.6,    [], 0,  0,   0,     5e5,    1e-2,       1e-6,   1};

    % Euler Resting --> Aligned Mixed12
% TestBatches(1).set =  {'ER2AM',    6.7,    [],  0,  1e-3,   .1,     1.5e6,    1e-2,       1e-6,   1};

    % Euler Resting --> Lagrange Re2sting
% TestBatches(1).set =  {'ER2LR',    1.8,    [],  0,  1e-3,   .1,     5e5,    1e-2,       1e-6,   1};
TestBatches(1).set =  {'ER2LR',    3.5,    [],  0,  0,   0,     5e5,    1e-2,       1e-6,   1};

    % Lagrange Resting --> Euler Resting
% TestBatches(1).set =  {'LR2ER',     5.2,    [],  0,  1e-3,   0,     50000,    1e-2,       1e-6,   1};
% TestBatches(1).set =  {'LR2ER',     1,    [],  0,  0,   0,     1e5,    1e-2,       1e-6,   1};

% ////////////////////////////////////////////////////////////

for i = 1:length(TestBatches)
    clear Body BodyOut tvec filename fullFilename
    Dura = TestBatches(i).set{7}; initStep = TestBatches(i).set{8};
    tol = TestBatches(i).set{9}; varstep = TestBatches(i).set{10};
    
    %Setup Filename
    datetime=datestr(now);
    datetime=strrep(datetime,':','_'); %Replace colon with underscore
    datetime=strrep(datetime,'-','_');%Replace minus sign with underscore
    datetime=strrep(datetime,' ','_');%Replace space with underscore
    filename = [TestBatches(i).set{1},'_',num2str(Dura/1000),'ks_tol',num2str(log10(tol)),'_',datetime];
    savelocation = '/Users/hartzell_lab/Documents/MATLAB/Research Project/2.0/Data Sets/Batch 2/';
        
    %Initialize Body Properties
    [Body, G] = testConfigsBatch(TestBatches(i).set);
    for a=1:3
        Body(a).contact(1,1:3)=0;
        Body(a).collision(1,1:3)=0;
    end

    tend = Dura/TU; 
    step = initStep/TU; %initial step size to try (or to use for fixed step propagation)
    fullFileName = fullfile(savelocation, filename);
    save(fullFileName)
    
    tinit = clock;
    [tvec, BodyOut] = mainPropagatorVerletBatch(tend,Body,G,varstep,step,tol,fullFileName); % Verlet
    % [tvec, BodyOut] = mainPropagator(tend,Body,G,varstep,step,tol); % RK variable timestep
    % [tvec, BodyOut] = SHAKE_prop(tend,Body,G,varstep,step,tol)
    tfinal = clock;
    telap = etime(tfinal,tinit);
    disp(['Elapsed time: ', num2str(telap),' seconds'])
    
    save(fullFileName,'-v7.3')
end


 
 