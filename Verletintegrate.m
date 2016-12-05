
function [nextState, flag, breakOut] = Verletintegrate(CurrState, h,G,N)

%%%%%% First check if bodies are separating %%%%%%
% breakOut = check_if_break(CurrState, G, N);
% orbitRates = calcOrbRates(CurrState, breakOut);
% CurrState.orbitRates = orbitRates;

[breakOut, ~] = check_if_break2(CurrState, h, G, N);

% % ^This step inadvertenly propagates the bodies if all contact points are broken, saving computation
% % time
% if any(any(CurrState.contact))
%     if ~any(~(breakOut.BC(:))) %checks if all contact points broke
%         nextState = nextStateTest;
%         return
%     end
% end
% breakOut = [];
[nextState, flag] = VVintegrate(CurrState, h, G, N, breakOut);
    
end