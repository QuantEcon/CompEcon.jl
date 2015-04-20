function demfin01

%% DEMFIN01
% Illustrate finites difference Hessian evaluation 

% Preliminary tasks
demosetup(mfilename)

figure
axis([-2 2 -2 2])
set(gca,'xtick',-2:1:2,'xticklabel','')
set(gca,'ytick',-2:1:2,'yticklabel','')
ftext(-1.0,-2,'x_1-h_1','center','top')
ftext( 0.0,-2,'x_1'    ,'center','top')
ftext( 1.0,-2,'x_1+h_1','center','top')
ftext(-2.3,-1,'x_2-h_2','center','middle')
ftext(-2.3, 0,'x_2'    ,'center','middle')
ftext(-2.3, 1,'x_2+h_2','center','middle')
ftext(-1, 1,'f^{ - +}','center','middle',18)
ftext( 0, 1,'f^{ 0 +}','center','middle',18)
ftext( 1, 1,'f^{ + +}','center','middle',18)
ftext(-1,-1,'f^{ - -}','center','middle',18)
ftext( 0,-1,'f^{ 0 -}','center','middle',18)
ftext( 1,-1,'f^{ + -}','center','middle',18)
ftext(-1, 0,'f^{ - 0}','center','middle',18)
ftext( 0, 0,'f^{ 0 0}','center','middle',18)
ftext( 1, 0,'f^{ + 0}','center','middle',18)
title('Evaluation Points for Finite Difference Hessians','fontsize',16)

% Save Plots as EPS Files
printfigures(mfilename,1)