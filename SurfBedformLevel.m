
function [] = SurfBedformLevel(X,Y,Z,H,Lvl,PrintPath)

% This function prints each level of the bedforms discrimination

% INPUTS
%     1. X, Y, Z: Coordinates for all the points in the bedform river
%     2. H: coordinate Z for the corresponding level,
%     3. Lvl: 1,2 or 3 depending of the level from the discrimination 
%     4. PrintPath: Path to the folder to save the figures
%%

f = figure();

if (Lvl>1)
    colormap hsv;
    surf(X,Y,H,'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud','SpecularStrength',0);
else
    %        blue, yellow, white
    ct = [ 0 0 1 ; 1 1 0 ; 1 1 1 ];
    colormap (ct);
    colorSet = Z;
    [x,y] = size(Z);
    
    minh = min(min(H));
    maxh = max(max(H));
    
    for i=1:x
        for j=1:y
            if H(i,j) <= -0.1 
                colorSet(i,j) = (-0.1 - H(i,j))/(-0.1-minh) ;
            elseif H(i,j) >= 0.1
                colorSet(i,j) = (0.1 - H(i,j))/(0.1-maxh) ;
            else
                colorSet(i,j) = (H(i,j)+0.1)/(0.2) ;
            end
        end
    end
    surf(X,Y,Z,colorSet,'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud','SpecularStrength',0);
end

Title  = sprintf(['Plot - Level ' num2str(Lvl)]);
% Title
set(get(gca,'Title') ,'String',Title);
set(get(gca,'Title') ,'FontSize',20);
set(get(gca,'Title') ,'FontWeight','bold');

% Axes labels
set(get(gca,'XLabel'),'String','X(m)');
set(get(gca,'XLabel'),'FontSize',18);
set(gca,'XTick',[]);
set(get(gca,'YLabel'),'String','Y(m)');
set(get(gca,'YLabel'),'FontSize',18);
set(gca,'YTick',[]);
set(gca,'ZTick',[]);

axis off

daspect([5 5 1]);
axis tight;
view(3);
camlight left;
if Lvl==1 
    view(2);
end

cm = colorbar('location','westoutside');

switch Lvl
    case 1  
        labelLvl = 'h_{1,3}';
    case 2  
        labelLvl = 'h_{2,3}';
    case 3  
        labelLvl = 'h_{3,3}';
end

ylabel(cm,[labelLvl ' [m]'],'fontsize',14,'fontweight','bold');

if (Lvl>1)
    set(cm, 'YTick',     []);
    set(cm, 'YTickLabel',[]);
else
    XT = [ 1/3 ; 2/3];
    XL = {'-0.1' ; '0.1'};
    set(cm, 'YTick',     XT);
    set(cm, 'YTickLabel',XL);
end
mTextBox = uicontrol('style','text','Units','normalized','FontSize',8,...
                    'String','Plotted by Bedforms-ATM','BackgroundColor',[1 1 1],...
                    'Position',[0.7 0.05 0.25 0.03]);
                
set(gcf,'paperunits','centimeters');
set(gcf,'papersize',[28,19]);
set(gcf,'paperposition',[0.05,0.05,29*.95,20*.95]); 
print(gcf,'-dpdf','-r300',[PrintPath Title]);
close(f);
