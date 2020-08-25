function [] = plotRelaxation(GDconc,GDflag,T1_stats,axisInfo,Report_folder)

GDplot = GDconc(GDflag);
T1plot = T1_stats(GDflag,1)/1000;
T1error = T1_stats(GDflag,2)/1000;

figure('rend','painters','pos',[10 10 300 300]);
subplot(2,1,1);
errorbar(GDplot, T1plot, T1error, 'o-', 'LineWidth',2); hold on;
xlim(axisInfo.xlim);
xlabel('Gd concentration (mM)'); ylabel(axisInfo.ylabel);

subplot(2,1,2);
plot(GDplot, 1./T1plot, 'o', 'LineWidth',2); hold on;
[F,gof] = fit(GDplot, 1./T1plot, 'poly1');
plot(F);
legend off;
text(0,6,['Relaxivity ' num2str(F.p1,'%.2f') ' s^{-1}mM^{-1}'],...
    'FontSize',14);
text(0,4.5,['R^2 ' num2str(gof.rsquare,'%.2f')],...
    'FontSize',14);
xlim(axisInfo.xlim); ylim([0 7]);
xlabel('Gd concentration (mM)'); ylabel(axisInfo.ylabel2);

export_fig([Report_folder '/' axisInfo.imgFile],'-m2','-transparent');
close;

end