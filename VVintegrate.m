
function [nextState, flag] = VVintegrate(CurrState, h, G, N, breakOut)

%%%%%% First check if bodies are separating %%%%%%
% breakOut = check_if_break(CurrState, G, N);
% orbitRates = calcOrbRates(CurrState, breakOut);
% CurrState.orbitRates = orbitRates;

flag.set=0;
v_mid=zeros(N,3); w_mid=zeros(N,3);

% % 1) Get a(t) = f(x(t))
% k1=calcForces2(CurrState,G,N,h, breakOut);
% 
% nextState=CurrState; %For matters of easily setting I, contact, and other parameters
% 
% % 2) Get x(t+h) = x(t) + v(t)*h + 1/2*a*h^2
% for i=1:N
%     nextState.r(i,:) = CurrState.r(i,:) + CurrState.v(i,:)*h + 1/2*k1.v(i,:)*h^2;
%     nextState.psi(i,1) = rem(CurrState.psi(i,1) + CurrState.w(i,3)*h + 1/2*k1.w(i,3)*h^2,2*pi);
% end
% 
% % 3) Get midstep velocity v(t+h/2)
% for i=1:N
%     v_mid(i,:) = CurrState.v(i,:) + h/2*k1.v(i,:);
%     w_mid(i,:) = CurrState.w(i,:) + h/2*k1.w(i,:);
% end
% 
% % 4) Now get a(t+h) = f(x(t+h), v(t+h/2))
% midState=nextState;
% midState.v=v_mid; midState.w=w_mid;
% k2=calcForces2(midState,G,N,h, breakOut);
% 
% % 5) Finally get v(t+h) = v(t+h/2) + 1/2*(a(t))*h
% for i=1:N
%    nextState.v(i,1:3) = v_mid(i,:) + h/2*k2.v(i,:);
%    nextState.w(i,1:3) = w_mid(i,:) + h/2*k2.w(i,:);
% end

% 1) Get a(t) = f(x(t))
k1=calcForces2(CurrState,G,N,h, breakOut);

nextState=CurrState; %For matters of easily setting I, contact, and other parameters

% 2) Get x(t+h) = x(t) + v(t)*h + 1/2*a*h^2
for i=1:N
    nextState.r(i,:) = CurrState.r(i,:) + CurrState.v(i,:)*h + 1/2*k1.v(i,:)*h^2;
    nextState.psi(i,1) = rem(CurrState.psi(i,1) + CurrState.w(i,3)*h + 1/2*k1.w(i,3)*h^2,2*pi);
end

% 3) Get midstep velocity v(t+h/2)
for i=1:N
    v_mid(i,:) = CurrState.v(i,:) + h/2*k1.v(i,:);
    w_mid(i,:) = CurrState.w(i,:) + h/2*k1.w(i,:);
end

% 4) Now get a(t+h) = f(x(t+h), v(t))
% midState=nextState;
% midState.v=v_mid; midState.w=w_mid;
k2=calcForces2(nextState,G,N,h, breakOut);

% 5) Finally get v(t+h) = v(t+h/2) + 1/2*(a(t))*h
for i=1:N
   nextState.v(i,1:3) = v_mid(i,:) + h/2*k2.v(i,:);
   nextState.w(i,1:3) = w_mid(i,:) + h/2*k2.w(i,:);
end

    
end