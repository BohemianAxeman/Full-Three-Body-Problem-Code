
function fixedState = fixCollision(CurrState, nextState, N)
global MU DU TU
%Organize body positions and velocities into easy to work with vector array
    r_vec(1:N,1:3)=CurrState.r(1:N,1:3);
    v_vec(1:N,1:3)=CurrState.v(1:N,1:3);
    w_vec(1:N,1:3)=CurrState.w(1:N,1:3);
    collision(1:N,1:N)=nextState.collision(1:N,1:N);
    m(1:N)=CurrState.m(1:N);
    R(1:N)=CurrState.R(1:N);
    Iinv = zeros(N,N,N);
    
    for j=1:N
        Iinv(1:N,1:N,j)=inv(diag(CurrState.Ivec(j,1:3)));
    end
e=0; %coefficient of elasticity

%First, move colliding bodies so they are just touching
fixedBodyPos=fixContactPos(CurrState, nextState, N);
r_vec(1:N,1:3)=fixedBodyPos.r(1:N,1:3);
contact=fixedBodyPos.contact;

for k=1:N-1
    for j=k+1:N
        dr=r_vec(j,1:3)-r_vec(k,1:3);
        n=-dr/norm(dr); %Normal vector, points from b(j) to a(k)
        vk=v_vec(k,:); %+cross(w_vec(k,:),dr/2); %Velocity of cotact point on body k
        vj=v_vec(j,:); %+cross(w_vec(j,:),-dr/2); %Velocity of contact point on body j
        vrel=vj-vk; %Relative velocity
        if (dot(vrel,n)>10^-6 * (TU/DU) && contact(k,j)) ||  collision(k,j) %% vrel > .1mm/s in direction of bodies
            collision(k,j)=1; collision(j,k)=1;
        end
    end
end

%Now do sequential 2-Body collision handling until some tolerance is met
while any(any(collision))
    for k=1:N-1
        for j=k+1:N
            
            if collision(k,j) %Relative velocity > 1mm/s, Collision has occured
                %Body k = a, Body j = b
                dr=r_vec(j,1:3)-r_vec(k,1:3);
                rap=dr/2; rbp=-dr/2; %Position of contact point relative to body a & b
                n=-dr/norm(dr); %Normal vector, points from b to a
                vap1=v_vec(k,:) ;%+ cross(w_vec(k,:),rap);
                vbp1=v_vec(j,:) ;%+ cross(w_vec(j,:),rbp);
                vab1=vap1-vbp1; %Relative velocity at collision point p
                term1=cross(Iinv(:,:,k)*cross(rap,n)', rap);
                term2=cross(Iinv(:,:,j)*cross(rbp,n)', rbp);
                jImp=-(1+e)*dot(vab1,n) / ( 1/m(k)+1/m(j) + dot(n,term1) + dot(n,term2)); %Calculate impulse factor j
                
                % Update translational and rotational velocities of bodies post collision
                    v_vec(k,:)=v_vec(k,:) + jImp/m(k)*n;
                    v_vec(j,:)=v_vec(j,:) - jImp/m(k)*n;
                % These below are the true equations for angular velocity post-collision
%                     w_vec(k,:)=w_vec(k,:) + (Iinv(:,:,k)*cross(rap,jImp*n)')';
%                     w_vec(j,:)=w_vec(j,:) - (Iinv(:,:,j)*cross(rbp,jImp*n)')';
                % For now, I use a simple approach to better approximate the rotation rates
                % post-collision assuming dynamic friction has caused
                % the two bodies to match angular velocities
                    %wf=(Iinv(3,3,k)*w_vec(k,:) - Iinv(3,3,j)*w_vec(j,:)) / (Iinv(3,3,k)+Iinv(3,3,j));
                    %w_vec(k,:)=wf; w_vec(j,:)=-wf;
                
                collision(k,j)=0; collision(j,k)=0;
            end %End of IF: collision phyics 
        end
    end %End of colliding bodies
    
    %Check for any persisting collisions
    for a=1:N-1 
        for b=a+1:N
            dr=r_vec(b,1:3)-r_vec(a,1:3); dv=v_vec(b,1:3)-v_vec(a,1:3);
            if contact(a,b) && dot(dv,-dr/norm(dr))>10^-12 * (TU/DU) %norm(dv)>10^-4 && dot(dv,-dr)>0
                collision(a,b)=1; collision(a,b)=1;
            end
        end
    end
end %End of WHILE

fixedState.r(1:N,1:3)=r_vec(:,:);
fixedState.v(1:N,1:3)=v_vec(:,:);
fixedState.w(1:N,1:3)=w_vec(:,:);
fixedState.contact=contact;
fixedState.m=CurrState.m;
fixedState.Ivec=CurrState.Ivec;

end