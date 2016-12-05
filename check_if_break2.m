% A function that propagates the current state of the system and sees whether the bodies in contact
% separate. If separation occurs, do not apply can spring/contact forces. Else if separation does
% not occur (bodies penetrate or hold surface contact) then apply spring forces.
function [breakOut, nextStateTest] = check_if_break2(CurrState, h, G, N)
    
    if any(any(CurrState.contact))
        testState = CurrState;
        testState.contact = zeros(N,N); %Assume no contact forces are needed
        [nextStateTest, ~] = VVintegrate(testState, h, G, N, []);
        [~, flag] = checkContact(CurrState, nextStateTest, N);
        nc = 0;
        for k=1:N-1
            for j=k+1:N
                if CurrState.contact(k,j)
                    nc = nc + 1;
                    breakOut.ContactPoints(nc).A = k;
                    breakOut.ContactPoints(nc).B = j;

                    if flag.brokeContact(k,j)
                        breakOut.BC(nc) = 1;
                    else
                        breakOut.BC(nc) = 0;
                    end
                end
            end
        end

        breakOut.nc = nc;
    else
        breakOut = [];
        nextStateTest = [];
    end
end