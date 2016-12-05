
% Plots the positions of bodies in an animation
% ** MUST RUN 'dataAnalysis' beforehand
clear t_slct
co = [    0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];
global DU TU MU
fig=figure;
hold all
xlabel('X position (km)')
ylabel('Y position (km)')
zlabel('Z position (km)')
% legend('Body 1','Body 2','Body 3')
% title('Body Positions')
axis equal

axis([-15 15 -15 15])
mTextBox = uicontrol('style','text');
set(mTextBox,'Units','characters')
set(fig, 'Position', [100, 100, 1000, 800]);
% set(mTextBox,'Position',[1 30 30 2])
mTextBox.FontSize=14;
grid on
% Divide simulation duration into equal-ish sections
t_slct=1; Tani=100; %animation time step
for i=2:length(t1)
    if t1(i)-t1(t_slct(end))>=Tani
        t_slct(end+1)=i;
    elseif i==length(t1)
        t_slct(end+1)=i; %Makes sure to include final time
    end
end

kwait=waitforbuttonpress;
% set(mTextBox,'Position',[1 41 30 2])
set(mTextBox,'Position',[131 48 30 2])
for i=1:length(t_slct)
    k=t_slct(i);
    h(1)=viscircles([r(k,1,1),r(k,2,1)] , radius, 'EdgeColor',co(1,:));
    harrow(1)=plot([r(k,1,1); r(k,1,1)+radius*cos(psi(k,1,1))],[r(k,2,1); r(k,2,1)+radius*sin(psi(k,1,1))],'Color',co(1,:));
    h(2)=viscircles([r(k,1,2),r(k,2,2)] , radius, 'EdgeColor',co(2,:));
    harrow(2)=plot([r(k,1,2); r(k,1,2)+radius*cos(psi(k,1,2))],[r(k,2,2); r(k,2,2)+radius*sin(psi(k,1,2))],'Color',co(2,:));
    h(3)=viscircles([r(k,1,3),r(k,2,3)] , radius, 'EdgeColor',co(3,:)); 
    harrow(3)=plot([r(k,1,3); r(k,1,3)+radius*cos(psi(k,1,3))],[r(k,2,3); r(k,2,3)+radius*sin(psi(k,1,3))],'Color',co(3,:));
    h(4) = viscircles([Rcm(k,1), Rcm(k,2)], radius/8, 'EdgeColor','k');
    harrow(4) = plot([Rcm(k,1)-radius/8; Rcm(k,1)+radius/8],[Rcm(k,2); Rcm(k,2)],'k');
    harrow(5) = plot([Rcm(k,1); Rcm(k,1)],[Rcm(k,2)-radius/8; Rcm(k,2)+radius/8],'k');
    set(mTextBox,'String',['Time = ' ,num2str(t1(k)),' s'])
    if i==1
        legend(h,'Body 1','Body 2','Body 3','Barycenter','Location','northeastoutside')
    end
    pause(.01)
    delete(h(:))
    delete(harrow(:))
end