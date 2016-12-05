function fixedBody = fixContactPos(CurrState, nextState, N)

fixedBody = CurrState;

%Organize Body Positions
r_vec(1:N,1:3)=CurrState.r(1:N,1:3);
v_vec(1:N,1:3)=CurrState.v(1:N,1:3);
R(1:N)=CurrState.R(1:N);
contact(1:N,1:N)=nextState.contact(1:N,1:N);
cntC=0;
for k=1:(N-1)
    for j=(k+1):N
        dr=r_vec(j,1:3)-r_vec(k,1:3);   
        if contact(k,j) || (norm(dr)-2*R(k))<=0 %Bodies have penetrated
            fixedBody.contact(k,j)=1; fixedBody.contact(j,k)=1;
            cntC=cntC+1;
        else %if norm(dr) > R(j)+R(k)
            fixedBody.contact(k,j)=0; fixedBody.contact(j,k)=0;
        end
    end
end
%Calculate Initial Center of Mass
% Rcmi=(r_vec(1,:)+r_vec(2,:)+r_vec(3,:))/3;
K=1;

if cntC == 1 % Only two bodies in contact
    for k=1:(N-1)
        for j=(k+1):N
            if fixedBody.contact(k,j)
            %Move both bodies equal amounts
            rrel=r_vec(j,:)-r_vec(k,:);
            r_vec(k,:)=r_vec(k,:)-K*(2*R(k)-norm(rrel))/2 * rrel/norm(rrel);
            r_vec(j,:)=r_vec(j,:)+K*(2*R(k)-norm(rrel))/2 * rrel/norm(rrel);
            end
        end
    end
elseif cntC == 2 % One body in contact with other two
    [~,a]=max(sum(fixedBody.contact,2)); %Central body
    temp=[1:N]; temp(a)=[];
    b=temp(1); c=temp(2);
    dr1=r_vec(b,:)-r_vec(a,:); dr2=r_vec(c,:)-r_vec(a,:);
    r_vec(b,:)=r_vec(b,:)+K*(2*R(k)-norm(dr1)) * dr1/norm(dr1);
    r_vec(c,:)=r_vec(c,:)+K*(2*R(k)-norm(dr2)) * dr2/norm(dr2);         
elseif cntC == 3 % All three bodies in contact
    a=1; %Choose body 1 as central body
    temp=[1:N]; temp(a)=[];
    b=temp(1); c=temp(2);
    dr1=r_vec(b,:)-r_vec(a,:); dr2=r_vec(c,:)-r_vec(a,:);
    r_vec(b,:)=r_vec(b,:)+K*(2*R(k)-norm(dr1)) * dr1/norm(dr1);
    r_vec(c,:)=r_vec(c,:)+K*(2*R(k)-norm(dr2)) * dr2/norm(dr2); 
end

% Rcmf=(r_vec(1,:)+r_vec(2,:)+r_vec(3,:))/3;
fixedBody.r(1:N,1:3)=r_vec(1:N,1:3);
end