function [] = plotCEST(pH,GLYconc,GLYflag,GDconc,GDflag,MTRasym_stats,...
    GLYarray,axisInfo,Report_folder)

GLYplot = GLYconc(GLYflag);
PHplot = pH(GLYflag);
MTRplot = MTRasym_stats(GLYflag,1);
MTRerror = MTRasym_stats(GLYflag,2);

figure('rend','painters','pos',[10 10 300 300]);
subplot(2,1,1);
for GLY = GLYarray
    errorbar(PHplot(GLYplot == GLY), ...
        MTRplot(GLYplot == GLY), MTRerror(GLYplot == GLY), ...
        'o-', 'LineWidth',2); hold on;
end
xlim(axisInfo.xlim);
xlabel('pH'); ylabel('MTRasym (%)');
legend(axisInfo.legend,'box','off');

% T2 dependency
GDflag = GDflag > 0;

GDplot = GDconc(GDflag);
MTRplot = MTRasym_stats(GDflag,1);
MTRerror = MTRasym_stats(GDflag,2);

subplot(2,1,2);
errorbar(GDplot, MTRplot, MTRerror, 'o-', 'LineWidth',2); hold on;
xlim(axisInfo.xlimGD);
xlabel('Gd concentration (mM)'); ylabel('MTRasym (%)');

export_fig([Report_folder '/MTRasymplot.png'],'-m2','-transparent');
close;

end