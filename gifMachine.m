clc
close all
clear all
% set(0,'DefaultFigureVisible','off')
nbr_cases = mainLoop(1);
f(1) = getframe(gcf);
close gcf

for i = 2:nbr_cases
    mainLoop(i);
    f(i) = getframe(gcf);
    close gcf
end
% set(0,'DefaultFigureVisible','on')

filename = 'testAnimated.gif'; % Specify the output file name
for i = 1:nbr_cases
    [A,map] = rgb2ind(f(i).cdata,256);
    if i == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',.25);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',.25);
    end
end

