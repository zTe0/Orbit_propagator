%% PLOTS
% folderpath = pwd;
close all

Sat_to_plot = str2num(app.SatellitetoplotninputnEditField.Value);
Walker_const = app.SwitchtoWalkerorbitCheckBox.Value; % e.g. GNSS
DATA = app.DATA;
DATA_SAT = app.DATA_SAT;
SampleTime = app.SampleTimesecondsEditField.Value;
SAT_search = app.QuerySatellitesCheckBox.Value;
ECEF_mode = app.ECEFCheckBox.Value;

ECI__ = app.ECI;
if ECEF_mode == 0
    LLA__ = app.LLA;
end

for j = 1:length(ECI__)
    if Sat_to_plot(1) > size(ECI__{j},1)
        Sat_to_plot(1) = size(ECI__{j},1);
    end
end

if Walker_const == false && SAT_search == false
    app.PlottedsatnameTextArea.Value = DATA.Name(Sat_to_plot(1));
elseif Walker_const == false && SAT_search == true
    app.PlottedsatnameTextArea.Value = DATA_SAT.Name(Sat_to_plot(1));
end

if length(ECI__) == 1
    Sat_to_plot(2) = 1;
end

planisphere(app)

for k = 1:2
    if k == 2
        figure(1)
        axes = gca;
    else
        axes = app.UIAxes1;
    end
if ECEF_mode == 0
    RF = 'ECI J2000';
    x = squeeze(LLA__{Sat_to_plot(2)}(Sat_to_plot(1),2,:))*180/pi;
    y = squeeze(LLA__{Sat_to_plot(2)}(Sat_to_plot(1),1,:))*180/pi;
    scatter(axes,x,y,15,'.r'), grid on, hold on
    xlim(axes,[-180 180]), ylim(axes,[-90 90])
    xticks(axes,-180:30:180), yticks(axes,-90:30:90)
    plot(axes,x(1),y(1),'og','markersize',8,'linewidth',3);
    plot(axes,x(end),y(end),'sqg','markersize',8,'linewidth',2.5);
    text(axes,x(1)+5,y(1)-5,'start','FontSize',16,'Color','g','FontWeight','bold')
    text(axes,x(end)-37,y(end)-5,'end','FontSize',16,'Color','g','FontWeight','bold')
    if Walker_const == false && SAT_search == false
        title(axes,sprintf('Ground Track (%s)',DATA.Name(Sat_to_plot(1))),'FontSize',18)
    elseif Walker_const == false && SAT_search == true
        title(axes,sprintf('Ground Track (%s)',DATA_SAT.Name(Sat_to_plot(1))),'FontSize',18)
    elseif Walker_const == true
        title(axes,'Ground Track (Walker-Delta)','FontSize',18)
    end
else
    RF = 'ECEF IAU2000';
    title(axes,'Ground Track unavailable','FontSize',18)
end
end

cd save/img
if ECEF_mode == 0
    savefig('GroundTrack.fig')
    % set(gcf, 'Position', get(0, 'Screensize'));
    % saveas(gcf,'GroundTrack.png')
    exportgraphics(gca,'GroundTrack.png')
    frm = "ECI";
else
    frm = "ECEF";
end

figure
if Walker_const == false
    if SAT_search == false
        T = seconds(minutes(DATA.T(Sat_to_plot(1))));
    else
        T = seconds(minutes(DATA_SAT.T(Sat_to_plot(1))));
    end
else
    T = app.T;
end
Earth3d(ECI__{Sat_to_plot(2)}(Sat_to_plot(1),:,:)*1e-3,T(Sat_to_plot(2)),SampleTime,ECEF_mode,app)
for k = 1:2
    if k == 2
        figure(2)
        axes = gca;
    else
        axes = app.UIAxes2;
    end
if Walker_const == false && SAT_search == false
    title(axes,sprintf('%s (%s)',RF,DATA.Name(Sat_to_plot(1))),'FontSize',18)
elseif Walker_const == false && SAT_search == true
    title(axes,sprintf('%s (%s)',RF,DATA_SAT.Name(Sat_to_plot(1))),'FontSize',18)
elseif Walker_const == true
    title(axes,sprintf('%s (Walker-Delta)',RF),'FontSize',18)
end
end

savefig(sprintf('%s.fig',frm))
% set(gcf, 'Position', get(0, 'Screensize'));
% saveas(gcf,sprintf('%s.png',frm))
exportgraphics(gca,sprintf('%s.png',frm))
cd ../..

% close figures, plots have been saved and are available in app only
close(1)
close(2)

%% Functions

function planisphere(app)

% Represents a map of the surface of Earth.

I = imread("EarthTex.jpg");
I = flip(I,1);

for k = 1:2
    if k == 1
        figure
        axes = gca;
    else
        axes = app.UIAxes1;
        hold(axes,'off')
    end
image(axes,[-180 180],[-90 90],I);
grid(axes,'on')
xticks(axes,-180:30:180)
yticks(axes,-90:30:90)
set(axes,'YDir','normal')
hold(axes,'on')
axis(axes,'equal','tight')
xlabel(axes,'Longitude [deg]','FontSize',15)
ylabel(axes,'Latitude [deg]','FontSize',15)
title(axes,'Ground Track','FontSize',18)
end
end

%--------------------------------------------------------------------------

function Earth3d(ECI,T,SampleTime,ECEF_mode,app)
xi = squeeze(ECI(:,1,1)); yi = squeeze(ECI(:,2,1)); zi = squeeze(ECI(:,3,1));
xf = squeeze(ECI(:,1,end)) ; yf = squeeze(ECI(:,2,end)); zf = squeeze(ECI(:,3,end));
if ECEF_mode == 1
    ECI = permute(ECI,[3,2,1]);
else
    ECI = permute(ECI(:,:,1:min([ceil(T/SampleTime)+1,size(ECI,3)])),[3,2,1]);
end
for k = 1:2
    if k == 1
        figure(2)
        axes = gca;
    else
        axes = app.UIAxes2;
        hold(axes,'off')
    end
% [X,Y,Z] = ellipsoid(0,0,0,6378.1366,6378.1366,6356.7519,50); % oblate
[X,Y,Z] = ellipsoid(0,0,0,6371,6371,6371,50); % spheric at arithmetic mean radius
I = imread("EarthTex.jpg");
props.FaceColor = 'texture';
props.EdgeColor = 'none';
props.CData = I;
surf(axes,X,Y,-Z,props), hold(axes,'on'), grid(axes,'on'), axis(axes,'equal'),% box on
plot3(axes,ECI(:,1),ECI(:,2),ECI(:,3),'LineWidth',1)
plot3(axes,xi,yi,zi,'og','MarkerSize',8,'linewidth',3)
plot3(axes,xf,yf,zf,'sqg','MarkerSize',8,'linewidth',2.5)
title(axes,'ECI J2000 (2BP)','FontSize',18)
xlabel(axes,'$x$[km]','Interpreter','latex','FontSize',18)
ylabel(axes,'$y$[km]','Interpreter','latex','FontSize',18)
zlabel(axes,'$z$[km]','Interpreter','latex','FontSize',18)
text(axes,xi-5,yi-5,zi-1000,'start','FontSize',16,'Color','g','FontWeight','bold')
text(axes,xf-5,yf-5,zf-1000,'end','FontSize',16,'Color','g','FontWeight','bold')
end
end