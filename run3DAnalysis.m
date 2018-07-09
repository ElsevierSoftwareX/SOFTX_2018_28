function totTime = run3DAnalysis (WLx, WLy, Xo, Yo, DLx, DLy, DataX, DataY,...
    DataZ, timetest, saveFigure, figureformat, outputPath)
    tic;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %CODE FOR ANALYSIS OF THE THREEDIMENSIONALITY INDEX Tb FOR THE
    %PARANÁ RIVER DATASET
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               %FRANCISCO NÚÑEZ-GONZÁLEZ 201604
    %==============================================
    
    %DataX and DataY starts in 0
    
    % Index of the LIMITS OF THE RANDOM FIELD TO ANALYZE
    Xmin = find(DataX(:,1) <= Xo , 1, 'last' );  %minimum CELL along x
    Xmax = find(DataX(:,1) >= Xo + DLx, 1) - 1;  %maximum CELL ALONG x
    Ymin = find(DataY(1,:) <= Yo , 1, 'last' );  %maximum CELL ALONG y
    Ymax = find(DataY(1,:) >= Yo + DLy, 1) - 1;  %maximum CELL ALONG y
    % SPECIFICATION OF THE LENGTH SCALE TO ANALYZE (in number of PAIR cells)
    Lx = find(DataX(:,1) <= WLx , 1, 'last' ); %Lx must be < (Xmax-Xmin) % RG: must be difined from the wavelet analysis
    Ly = find(DataY(1,:) <= WLy , 1, 'last' ); %Lx;  %Ly must be < (Ymax-Ymin) % RG: must be difined from the wavelet analysis
    
    %==============================================
    %_________________________ 
    %ANALYSIS FOR EACH NODE
    % RG/May16-2016: Here, it is needed to get the minimum location to avoid data
    % normalization
    %First and last node to analyze along x
    Xini = Xmin+ceil(Lx/2);
    Xend = Xmax-ceil(Lx/2);
    %First and last node to analyze along y
    Yini = Ymin+ceil(Ly/2);
    Yend = Ymax-ceil(Ly/2);
    %_________________________   
    %Loop for each cell
    it = 0; itx = 0; % RG sets the global counts for iterations
    % totIts = (1+Xend-Xini)*(1+Yend-Yini); % RG defines the total number of iterations

    % Inserted by RG ******************************************
    lenItX = Xend - Xini + 1;
    lenItY = Yend - Yini + 1;
    xTbm = zeros(lenItX,lenItY);
    yTbm = zeros(lenItX,lenItY);
    InTbm = zeros(lenItX,lenItY);            
    % Inserted by RG ******************************************

    % RG: these loops perform the analysis in stepSize = 1
    if ~timetest
        wb = waitbar(0,'1','Name','Progress...');
    end
    for xRel = Xini:Xend
    itx = itx+1;
    ity = 0;
       for yRel = Yini:Yend % RG relative Y coordinate
          if ~timetest
            waitbar(it/((1+Xend-Xini)*(1+Yend-Yini)),wb,sprintf('Progress in %.0f %%',round(100*it/((1+Xend-Xini)*(1+Yend-Yini)),0)));
          end
          it = it+1;
          ity = ity+1;
          %Select the frame for the analyzed node
          xi = xRel-ceil(Lx/2); xf = xRel+ceil(Lx/2); %Initial and final x cordinate
          yi = yRel-ceil(Ly/2); yf = yRel+ceil(Ly/2); %Initial and final y coordinate (transept number)       
          % Inserted by RG *******************************************

          xWindow = DataX ( xi:xf, yi:yf )';
          zWindow = DataZ ( xi:xf, yi:yf )';
          %STATISTICAL ANALYISIS OF SELECTED DATA FRAME
          %Size of the frame
          [lonYWindow,lonXWindow] = size(xWindow); %along x 
          % = cy;          %along y % RG: Change this variable name ??? <========
          zmeany = mean(zWindow(1:lonYWindow,1:lonXWindow));
          zmeanx = mean(zWindow(1:lonYWindow,1:lonXWindow),2)';
          zmeanxy = mean2(zWindow(1:lonYWindow,1:lonXWindow));
          STDx = std(zWindow(1:lonYWindow,1:lonXWindow),1,2)';
          STDxy = std2(zWindow(1:lonYWindow,1:lonXWindow));


          % Added by RG *****************************************                    
          zTemp = zeros(lonYWindow,1);
          for j=1:lonYWindow
            zTemp(j) = sum((zWindow(j,1:lonXWindow)-zmeanx(j))./...
                (STDx(j)).*((zmeany(1:lonXWindow)-zmeanxy)/STDxy));  
          end 
          R = sum(zTemp)/(lonYWindow*lonXWindow(1));            
          % Added by RG *****************************************
          Tb = 1-R^2; %THREE-DIMENSIONALITY INDEX 
          xTbm(itx,ity) = xRel;
          yTbm(itx,ity) = yRel;
          InTbm(itx,ity) = Tb;
          if timetest
            totTime = toc*(1+Xend-Xini)*(1+Yend-Yini);
            return
          end
       end
    end
    close (wb);
        titleFont = 15;
        axisFont = 14;
        labelFont = 15;
        figure();
        axes1 = axes('fontsize',axisFont);
        colorbar('location','eastoutside')
        hold on
        contourf(xTbm,yTbm,InTbm);
        Serie = ['Lx=' num2str(round(WLx,0)) ' m, Ly=' num2str(round(WLy,0)) ' m']; 
        title(axes1,Serie,'fontsize',titleFont,'FontWeight','bold');
        xlabel(axes1,'x (m)','fontsize',labelFont);
        ylabel(axes1,'y (m)','fontsize',labelFont);
        zlabel('Tb');
        box(axes1,'on'); 
        axis equal
        uicontrol('style','text','Units','normalized','FontSize',8,...
                    'String','Plotted by Bedforms-ATM','BackgroundColor',[1 1 1],...
                    'Position',[0.7 0.05 0.25 0.03]);
        if saveFigure
                       
            base_figureName = [outputPath ' - Lx ' num2str(round(WLx,0)) ' Ly ' num2str(round(WLy,0))];
            
            pos = 1;
            while( exist([base_figureName '_' num2str(pos) '.' figureformat ],'file'))
                pos = pos + 1;
            end
            figureName = [base_figureName '_' num2str(pos) ];
            
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
                    saveas(gcf,figureName,'fig');
            end
     
        end
end