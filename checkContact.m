function [BodyOut, flag] = checkContact(CurrState, nextState, N)

global MU DU TU
flag.set=0; flag.newContact=zeros(N,N); flag.newCollision=zeros(N,N); flag.brokeContact=zeros(N,N);
r_vec(1:N,1:3)=nextState.r(1:N,1:3);
v_vec(1:N,1:3)=nextState.v(1:N,1:3);
R(1:N)=nextState.R(1:N);
contact=CurrState.contact(:,:); %Previous contact condition
collision=zeros(N,N);

%Check contact/collision condition of nextState
for k=1:(N-1) %N-1
    for j=(k+1):N %N
        dr=r_vec(j,1:3)-r_vec(k,1:3);
        %Check if the two bodies are touching with some significant
        %relative velocity
        if norm(dr)<=(R(j)+R(k)) %Bodies are in contact
            if contact(k,j)==0 %New Contact formed
                contact(k,j)=1; contact(j,k)=1;
                flag.set=1;
                flag.newContact(k,j)=1; flag.newContact(j,k)=1;
            end
            
            n=-dr/norm(dr); %Normal vector, points from b(j) to a(k)
            vk=v_vec(k,:); %+cross(w_vec(k,:),dr/2); %Velocity of cotact point on body k
            vj=v_vec(j,:); %+cross(w_vec(j,:),-dr/2); %Velocity of contact point on body j
            vrel=vj-vk; %Relative velocity
            
            if dot(vrel,n)> 10^-6 * (TU/DU) %% vrel > 1mm/s * (TU/DU) in direction of bodies
                % account collision condition for Canonical units
                %'Bodies have collided!'
                collision(k,j)=1; collision(j,k)=1;
                flag.newCollision(k,j)=1; flag.newCollision(j,k)=1;
                flag.set=2;
            end 
        elseif contact(k,j) %Bodies not in contact but contact condition is set
            % Contact has broken
            flag.brokeContact(k,j)=1; flag.brokeContact(j,k)=1; 
            contact(k,j) = 0; contact(j,k) = 0;
            %flag.set=-1;
        end
        
    end
end
BodyOut.contact=contact; %This is the contact condition of nextState
BodyOut.collision=collision;

end