% Function that Plots all the figures
function plotFiguresWlt(motherWavelet, wltParameter, titleFigA,...
        signalDataX, signalData, wltPower, fourierPeriod, ...
        coneOfInfluence, wltSignif, globalSignif, globalWs, lag1,...
        saveFigure, figureformat, outputPath, plot_a, plot_b, plot_c,...
        signalName, signalUnit, sampleName, sampleUnit,...
        peaksVal, xPeakspos, yPeakspos)

    % This function prints each profile and its corresponding wavelet
    % spectrum and significance.

    % INPUTS
    %     Called from wavelet Analysiss
    % 
    %%
    
% Text size for the plot
titleFont = 15;
labelFont = 15;
axisFont = 14;

% Variables
signalVariance = std(signalData)^2;

% Range of the plots
plotXRange = [min(signalDataX), max(signalDataX)];

% Ticks of the Fourier Period
Yticks = 2.^(fix(log2(min(fourierPeriod))):fix(log2(max(fourierPeriod))));

% Plot Figure (a)
subplot(plot_a);

plot(signalDataX, signalData, 'b','linewidth',0.8);
ylabel(gca,[signalName ' (' signalUnit ')'],'fontsize',labelFont);
title (strcat('(a)', {' '}, titleFigA),'fontsize',titleFont,'FontWeight','bold');
set(gca, 'XLim', plotXRange(:));
set(gca, 'fontsize',axisFont);
grid on;
hold on;
plot(signalDataX,mean(signalData),'red--','linewidth',0.4);
hold off;

% Plot Figure (b)

subplot(plot_b);

titleFigB=strcat('(b)', {' '}, motherWavelet, {' '}, num2str(wltParameter),...
    {' '}, 'wavelet power spectrum');
Yticks = 2.^(fix(log2(min(fourierPeriod))):fix(log2(max(fourierPeriod))));
% 95% confidence contours
contour(signalDataX,log2(fourierPeriod),wltSignif,[-99,1],'k','linewidth',1);
xlabel(gca,[sampleName ' (' sampleUnit ')'],'fontsize',labelFont);
ylabel(gca,['\lambda (' signalUnit ')'],'fontsize',labelFont);
set(gca, 'XLim', plotXRange(:));
set(gca, 'YLim', log2([min(fourierPeriod),max(fourierPeriod)]),...
    'YDir','reverse', 'YTick',log2(Yticks(:)), 'YTickLabel',Yticks);
title(titleFigB,'fontsize',titleFont,'FontWeight','bold');
set(gca,'fontsize',axisFont);
hold on;

% Cone of influence
plot(signalDataX,log2(coneOfInfluence),'blue--','linewidth',2.0) 
%colorbar('location','eastoutside');
box on;
hold off;

% Plot Figure (c) % Originally it was Figure (e)
subplot(plot_c);

semilogx(globalWs,log2(fourierPeriod),'linewidth',1.1);
hold on;
grid on;
semilogx((globalSignif),log2(fourierPeriod),'red--');
xlabel(gca,['Ave. variance (' signalUnit '^2 )'],'fontsize',labelFont);
ylabel(gca,['\lambda (' signalUnit ')'],'fontsize',labelFont);
set(gca,'YLim',log2([min(fourierPeriod),max(fourierPeriod)]),...
    'YDir','reverse', 'YTick',log2(Yticks(:)), 'YTickLabel',Yticks);
set(gca,'fontsize',axisFont);
% Fourier Peaks
lengthWS = length(peaksVal);


for r=1:lengthWS
    text(xPeakspos(r), yPeakspos(r),num2str(peaksVal(r),'%6.1f'),'fontsize',13);
end

hold off;
title (strcat('(c) Global WT spectrum', {' '}, '\alpha','=',...
    num2str(lag1,2)), 'fontsize',titleFont,'FontWeight','bold');

% Save the Figure
if saveFigure
        figureName = [outputPath filesep strcat(strtrim(titleFigA),'_', motherWavelet, num2str(wltParameter))];
        set(gcf,'paperunits','centimeters');
        set(gcf,'papersize',[28,19]);
        set(gcf,'paperposition',[0.05,0.05,29*.95,20*.95]);  
        switch figureformat
            case 'PDF'
                print(gcf, '-dpdf',figureName)
            case 'jpg'
                print(gcf,'-djpeg',figureName)
            case 'tiff'
                print(gcf,'-dtiff',figureName)
            case 'fig'
                set(gcf,'Visible','on');
                saveas(gcf,figureName,'fig');
                set(gcf,'Visible','off');
        end
     
end;

