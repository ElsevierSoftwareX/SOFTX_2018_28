function [] = plotBFS( signalAbcise, h, h13, h23, h33, name, figureformat,...
                       outputpath, sampleName, sampleUnit  )
    
    % This function is called to build the figures of the bedforms
    % discrimination.

    % INPUTS
    %     1. signalAbcise: Horizontal coordinate of the profile,
    %     2. h: Horizontal coordinate of the profile,
    %     3. h13, h23, h33: Vertical coordinates of each level,
    %     4. name, figure format, outputpath, sampleName, sampleUnit: 
%            parameters used to print the figure
    % This values are outputs of the statisticsWlt function
    % 
    %%
    v = figure('Visible','off');
    titleFont = 15;
    labelFont = 15;
    axisFont = 12;
    if ( strcmp(figureformat,'fig') ) 
        figure();
    else
        figure('Visible','off');
    end
    XLimits =  [min(signalAbcise) max(signalAbcise)];
    XSize = length(signalAbcise);
    % Original profile
    subplot('position',[0.1 0.7 0.85 0.25]);
    plot(signalAbcise, h,'b');
    hold on;
    plot(signalAbcise,ones(XSize,1) * mean(h),'r--','linewidth',0.4);
    ylabel('h','fontsize',titleFont,'FontWeight','bold');
    y = get(gca,'YTick');
    set(gca,'YTickLabel', cellstr(num2str(y(:), '%.2f')) )
    set(gca, 'XTick', []);
    set(gca, 'XLim', XLimits);
    set(gca, 'fontsize',axisFont);
    axis tight;

    % 1st Level
    YRange  = 0.1*(max(h13)-min(h13));
    YLimits = [ min(h13)-YRange max(h13)+YRange ];
    subplot('position',[0.1 0.5 0.85 0.18]);
    plot(signalAbcise,h13,'b');
    hold on;
    plot(signalAbcise,ones(XSize,1) * mean(h13),'r--','linewidth',0.4);
    ylabel('h_{1,3}','fontsize',labelFont,'FontWeight','bold');
    set(gca, 'XTick', []);
    y = get(gca,'YTick');
    set(gca,'YTickLabel', cellstr(num2str(y(:), '%.2f')) )
    set(gca, 'XLim', XLimits);
    set(gca, 'YLim', YLimits);
    set(gca, 'fontsize',axisFont);
    axis tight;
    
    % 2nd Level
    YRange  = 0.1*(max(h23)-min(h23));
    YLimits = [ min(h23)-YRange max(h23)+YRange ];
    subplot('position',[0.1 0.3 0.85 0.18]);
    plot(signalAbcise,h23,'b');
    hold on;
    plot(signalAbcise,mean(h23),'r--','linewidth',0.4);
    ylabel('h_{2,3}','fontsize',labelFont,'FontWeight','bold');
    set(gca, 'XTick', []);
    y = get(gca,'YTick');
    set(gca,'YTickLabel', cellstr(num2str(y(:), '%.2f')) )
    set(gca, 'XLim', XLimits);
    set(gca, 'YLim', YLimits);
    set(gca,'fontsize',axisFont);
    axis tight;

    % 3rd Level
    YRange  = 0.1*(max(h33)-min(h33));
    YLimits = [ min(h33)-YRange max(h33)+YRange ];
    subplot('position',[0.1 0.1 0.85 0.18]);
    plot(signalAbcise, h33,'b','linewidth',0.4);
    xlabel([sampleName '(' sampleUnit ')'],'fontsize',labelFont,'FontWeight','bold');
    ylabel('h_{3,3}','fontsize',labelFont,'FontWeight','bold');
    set(gca,'fontsize',axisFont);
    y = get(gca,'YTick');
    set(gca,'YTickLabel', cellstr(num2str(y(:), '%.2f')) )
    set(gca, 'XLim', XLimits);
    set(gca, 'YLim', YLimits);
    axis tight;
    hold off;
    
    [~,~,~] = mkdir( [outputpath 'Figures' ] );
    
    figureName = [outputpath 'Figures' filesep name '-BFPSeparation'];
    if ( exist([figureName '.' figureformat],'file') )
        delete([figureName '.' figureformat]);
    end
    % Sets the paper size and plots the file
    set(gcf,'paperunits','centimeters')
    set(gcf,'papersize',[28,19]) % original [28,19] 
    set(gcf,'paperposition',[0.05,0.02,29*.95,20*.95]) % [0.05,0.05,29*.95,20*.95]
    switch figureformat
        case 'PDF'
            print(gcf,'-dpdf',figureName)
        case 'jpg'
            print(gcf,'-djpeg',figureName)
        case 'tiff'
            print(gcf,'-dtiff',figureName)
        case 'fig'
            saveas(gcf,figureName,'fig');
    end
    close(v);