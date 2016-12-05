
function [tvec, BodyOut] = mainPropagatorVerlet(tend,Body,G,varstep,step,tol)
% dbstop if error
BodyOut = Body;
N = length(Body);
i = 1; %Succesful propagations
iteration = 0; %Propagation attempts
it_repeat = 0;
tvec = 0;
h = step; %Initial step size
havg = 0;
maxI = 500;
global currTime

%%%%%%%%%%%%%%%%%%%% Initial Check for Contacts/Collisions %%%%%%%%%%%%%%%%%%%%
    resetState.collision=zeros(N,N); resetState.contact=zeros(N,N);
    for a=1:N
        CurrState.r(a,1:3)=BodyOut(a).r(1,:);
        CurrState.v(a,1:3)=BodyOut(a).v(1,:);
        CurrState.w(a,1:3)=BodyOut(a).w(1,:);
        CurrState.R(a)=BodyOut(a).R;
        CurrState.I(a)=BodyOut(a).I; CurrState.Ivec(a,1:3)=BodyOut(a).Ivec;
        CurrState.m(a)=BodyOut(a).mass;
        CurrState.semiaxis(a,1:3)=BodyOut(a).semiaxis(1:3);
    end %Transfer initial conditions to CurrState
    fixedBody=fixCollision(CurrState, resetState, N);
    for a=1:N %Update initial Body states
        BodyOut(a).contact(1,1:N)=fixedBody.contact(a,1:N);
        BodyOut(a).r(1,1:3)=fixedBody.r(a,1:3);
        BodyOut(a).v(1,1:3)=fixedBody.v(a,1:3);
    end
tic

if varstep == 0 %fixed step size, preallocate memory
    numSteps = floor(tend/step)+1;
    tvec = zeros(numSteps,1);
    for a=1:N
        BodyOut(a).r(2:numSteps,1:3)=zeros(numSteps-1,3);
        BodyOut(a).v(2:numSteps,1:3)=zeros(numSteps-1,3);
        BodyOut(a).psi(2:numSteps,1)=zeros(numSteps-1,1);
        BodyOut(a).w(2:numSteps,1:3)=zeros(numSteps-1,3);
        BodyOut(a).contact(2:numSteps,1:N)=zeros(numSteps-1,3);
    end
end
%%%%%%%%%%%%%%%%%%%% Main Propagation Loop %%%%%%%%%%%%%%%%%%%%    
while tvec(i) < tend
    havg = havg + h;
    currTime=tvec(i);
    if ~mod(iteration,maxI)
        elapt = toc;
        havg = havg/maxI;
        disp('[Iteration #], Percent Complete, Current Time: ')
        disp(strcat('[',num2str(iteration),'], ', num2str(tvec(i)/tend*100),'% ,', num2str(tvec(i)))) %Percent propagation completed
        disp(['Estimated time to Completion: ',num2str((tend-tvec(i))/h*(elapt/maxI)/60),' min'])
        tic
    end
    iteration=iteration+1;
    t_int(iteration)=currTime;

    if iteration>1
        if t_int(iteration)==t_int(iteration-1)
            it_repeat=it_repeat+1;
            if it_repeat>=100
                %timestep has been stuck for 100 attempts
                disp('Propagator stuck')
                % Output results obtained up to this point
                tvec=tvec';
                return
            end
        else it_repeat=0;
        end
    end %Checks if propagator gets stuck
        
    %Option to stop propagation exactly at input end time
    if (tvec(i)+h)>tend %Overstep past tend, reduce h
        h=tend-tvec(i); %ensure next time step propagates to tend
    end
    
    for a=1:N %Transfer Current iteration Body states to CurrState
        CurrState.r(a,1:3)=BodyOut(a).r(i,:);
        CurrState.v(a,1:3)=BodyOut(a).v(i,:);
        CurrState.psi(a,1)=BodyOut(a).psi(i,1);
        CurrState.w(a,1:3)=BodyOut(a).w(i,:);
        CurrState.contact(a,1:N)=BodyOut(a).contact(i,1:N);
        CurrState.R(a)=BodyOut(a).R;
        CurrState.I(a)=BodyOut(a).I; CurrState.Ivec(a,1:3)=BodyOut(a).Ivec;
        CurrState.m(a)=BodyOut(a).mass;
        CurrState.semiaxis(a,1:3)=BodyOut(a).semiaxis(1:3);
    end
    
    %%%%%%%%%% Propgagte CurrState by h/2 (for variable step only) %%%%%%%%%%
    if varstep == 1
%         [breakOut, ~] = check_if_break2(CurrState, h/2, G, N);
%         [nextState1, ~] = VVintegrate(CurrState, h/2, G, N, breakOut);
        [nextState1, ~] = Verletintegrate(CurrState, h/2, G, N);

        %%%%%%%%%% Check for new contact/collsion %%%%%%%%%%
        [~, flag1] = checkContact(CurrState, nextState1, N);

        if flag1.set ~= 0
            tempState.contact = flag1.newContact;
            tc = Time2Collision(CurrState,tempState,N); 
            if tc > step
                % Propagation time large, must propagate normally
                h = .95*tc; %undershoot collision time
                continue
            elseif tc<0
            %Collision/contact occured in past, go back and try again
                h=h/2;
                continue
            elseif tc < step && tc~=0
                %tc is either 0 or < step, propagate with simple Euler step
                flag3.set=-1;
                while flag3.set~=0
                    [nextState1, ~]=Verletintegrate(CurrState,tc,G,N); 
                    [~, flag3]=checkContact(CurrState, nextState1, N);
                    if flag3.set~=0
                        tc=tc*.98;
                    end
                end
    % Don't propagate if tc = 0
                i=i+1; tvec(i)=tvec(i-1)+tc;
                for a=1:N %Update BodyOut
                    BodyOut(a).r(i,1:3)=nextState1.r(a,1:3);
                    BodyOut(a).v(i,1:3)=nextState1.v(a,1:3);
                    BodyOut(a).psi(i,1)=nextState1.psi(a,1); 
                    BodyOut(a).w(i,1:3)=nextState1.w(a,1:3);  
                end
            end

            if flag1.set == 1 %New contact occured, update contact info
                resetState.contact=flag1.newContact;
                fixedBody=fixContactPos(nextState1, resetState, N);
                for a=1:N
                    BodyOut(a).contact(i,1:N)=fixedBody.contact(a,1:N);
                    BodyOut(a).r(i,1:3)=fixedBody.r(a,1:3);
                end

            elseif flag1.set == 2 %Collision occured, apply collision physics to Bodies
                resetState.collision=flag1.newCollision; resetState.contact=flag1.newContact;
                fixedBody=fixCollision(nextState1, resetState, N);
                for a=1:N
                    BodyOut(a).contact(i,1:N)=fixedBody.contact(a,1:N);
                    BodyOut(a).r(i,1:3)=fixedBody.r(a,1:3);
                    BodyOut(a).v(i,1:3)=fixedBody.v(a,1:3);
                end
            end
            h=step;
            continue
        end
   
        %%%%%%%%%% Propgagte nextState1 by h/2 (for variable step only) %%%%%%%%%%
%         [nextState2, ~] = VVintegrate(nextState1, h/2, G, N, breakOut);
        [nextState2, ~] = Verletintegrate(nextState1, h/2, G, N);

        %%%%%%%%%% Check for new contact/collsion %%%%%%%%%%
        [~, flag2] = checkContact(CurrState, nextState2, N);

        if flag2.set ~= 0
            % Need smaller timestep to trigger collision in first h/2 step
            h = h*.75;
            continue
        end
        
    end
    %%%%%%%%%% Propgagte CurrState by h %%%%%%%%%%
%         [breakOut, ~] = check_if_break2(CurrState, h, G, N);
%         [nextState, ~] = VVintegrate(CurrState, h, G, N, breakOut);
        [nextState, ~] = Verletintegrate(CurrState, h, G, N);

        %%%%%%%%%% Check for new contact/collsion %%%%%%%%%%
        [~, flag] = checkContact(CurrState, nextState, N);

        if flag.set ~= 0
            tempState.contact = flag.newContact;
            tc = Time2Collision(CurrState,tempState,N); 
            if tc > step
                % Propagation time large, must propagate normally
                h = .95*tc; %undershoot collision time
                continue
            elseif tc<0
            %Collision/contact occured in past, go back and try again
                h=h/2;
                continue
            else
                %tc is either 0 or < step, propagate with simple Euler step
                flag3.set=-1;
                while flag3.set~=0
                    [nextState, ~]=Verletintegrate(CurrState,tc,G,N); 
                    [~, flag3]=checkContact(CurrState, nextState, N);
                    if flag3.set~=0
                        tc=tc*.98;
                    end
                end

                i=i+1; tvec(i)=tvec(i-1)+tc;
                for a=1:N %Update BodyOut
                    BodyOut(a).r(i,1:3)=nextState.r(a,1:3);
                    BodyOut(a).v(i,1:3)=nextState.v(a,1:3);
                    BodyOut(a).psi(i,1)=nextState.psi(a,1); 
                    BodyOut(a).w(i,1:3)=nextState.w(a,1:3);  
                end
            end

            if flag.set == 1 %New contact occured, update contact info
                resetState.contact=flag.newContact;
                fixedBody=fixContactPos(nextState1, resetState, N);
                for a=1:N
                    BodyOut(a).contact(i,1:N)=fixedBody.contact(a,1:N);
                    BodyOut(a).r(i,1:3)=fixedBody.r(a,1:3);
                end

            elseif flag.set == 2 %Collision occured, apply collision physics to Bodies
                resetState.collision=flag.newCollision; resetState.contact=flag.newContact;
                fixedBody=fixCollision(nextState1, resetState, N);
                for a=1:N
                    BodyOut(a).contact(i,1:N)=fixedBody.contact(a,1:N);
                    BodyOut(a).r(i,1:3)=fixedBody.r(a,1:3);
                    BodyOut(a).v(i,1:3)=fixedBody.v(a,1:3);
                end
            end
            h=step;
            continue
        end
    
%     %%%%%%%%%% Compute next time step %%%%%%%%%%
    if varstep==0 %Constant timestep
        hnext=step;
    else
        maxError=zeros(N,1);
        for a=1:N
            rError=max(abs(nextState.r(a,1:3)-nextState2.r(a,1:3)));
            vError=max(abs(nextState.v(a,1:3)-nextState2.v(a,1:3)));
            psiError=max(abs(nextState.psi(a,1)-nextState2.psi(a,1))); 
            wError=max(abs(nextState.w(a,1:3)-nextState2.w(a,1:3))); 
            maxError(a)=max([rError, vError, psiError, wError]);
        end
        delta=max(maxError); %Pull out max error in difference
        if delta==0
            %Error is immaculate, double h (handles case where delta is 0)
            hnext=2*h;
        else    
            h0 = .9*h*abs(tol/delta)^.1; %Change to h based on max state variable error
            hnext = h0;
            %%% Compare Acceptable Erorr %%%
            if abs(delta) > tol %Difference in states >= allowable tolerance
                %Propagation failed. Recompute h and redo previous propagation
                h=hnext;
                continue %Continue without updating states
            end
        end
    end
        
    %%%%%%%%%% Prepare next state %%%%%%%%%%
    i = i+1; tvec(i) = tvec(i-1)+h;
    h=hnext;

    resetState.contact=zeros(N,N);
    fixedBody=fixContactPos(nextState, resetState, N);
    nextState.contact = fixedBody.contact;
    
    for a=1:N %Update BodyOut
        BodyOut(a).r(i,1:3)=nextState.r(a,1:3);
        BodyOut(a).v(i,1:3)=nextState.v(a,1:3);
        BodyOut(a).psi(i,1)=nextState.psi(a,1); 
        BodyOut(a).w(i,1:3)=nextState.w(a,1:3);  
        BodyOut(a).contact(i,1:N)=nextState.contact(a,1:N);
    end
    
    
    
end %End of While Loop

tvec=tvec';

end
