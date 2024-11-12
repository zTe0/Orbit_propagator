clc, %clear, close all
% folderpath          = app.DirectoryPathEditField.Value;
% folderpath = cd;
% cd(folderpath)

evalin('base','clear')
[~,folder_name] = fileparts(pwd);
folderpath = pwd;
app.folderpath = folderpath;
cd ../
addpath(genpath(folder_name))
cd (folder_name)
tic0 = tic;
tic


%% DATA
% ASSIGN APP DATA TO MATLAB VARIABLES
Walker_const        = app.SwitchtoWalkerorbitCheckBox.Value; % e.g. GNSS
Sim_Time            = str2num(app.SimulationTimedurationEditField.Value); %
SampleTime          = app.SampleTimesecondsEditField.Value; % [s]
Periodcycle_mode    = app.PeriodiccyclemodeCheckBox.Value;
    TTT             = str2num(app.PeriodEditField.Value); % period to cycle all walker orbits, gives warning if higher than automatic value
ECEF_mode           = app.ECEFCheckBox.Value;
AUTO_period_mod     = app.AutoCheckBox.Value;

if Walker_const == false
Const_name          = app.ConstellationNameEditField.Value;
Max_sat             = app.NumberofsatellitesperpageEditField.Value;  % n of Sats per Sat page, number or 'All'
    if ~isempty(str2num(Max_sat))
        Max_sat = str2num(Max_sat);
    end
htr                 = str2num(app.htrkmEditField.Value); % minimum altitude threshold
p                   = app.PageEditField.Value; % Sat page number
update_TLE          = app.updateTLEonlineCheckBox.Value;  % boolean, update TLE from website
SAT_search          = app.QuerySatellitesCheckBox.Value;
Sat_name            = app.SatelliteSearchEditField.Value; %

else
    i               = str2num(app.idegEditField.Value); % [deg] inclination
    TT_P            = str2num(app.TnofsatsEditField.Value); % total number of satellites per plane %%, must be multiple of P
    P               = str2num(app.PplanesEditField.Value);  % number of equally spaced geometric planes 
    F               = str2num(app.FspacingEditField.Value);  % spacing between satellites in adjacent planes
    a               = str2num(app.akmEditField.Value); % [km] semimajor axis
    e               = str2num(app.eEditField.Value); % eccentricity
    om              = str2num(app.degEditField.Value);  % [deg] ArgOfPeriapsis
    StartDate       = str2num(app.EpochEditField.Value); % Epoch
end
free_memory         = app.MaxarraysizeBytesEditField.Value;
% folderpath          = app.DirectoryPathEditField.Value;
savevarname         = app.LLAvariablenameEditField.Value;  % Large file will be read after being stored
Sat_to_plot         = str2num(app.SatellitetoplotninputnEditField.Value);

% Ground Station
GS.Lat  = str2num(app.LatitudegeodeticradEditField.Value);   % [rad]
GS.Lon  = str2num(app.LongituderadEditField.Value);          % [rad]
GS.h    = str2num(app.AltitudegeodetickmEditField.Value);    % [m]
GS.EL_min  = str2num(app.MinimumElevationradEditField.Value);   % [rad]
GS.Name = app.GroundStationNameEditField.Value;
coord_es = [GS.Lat GS.Lon GS.h];

%--------------------------------------------------------------------------
%% RETRIEVAL
if Walker_const == false

filename = sprintf('%s.tle',Const_name);
filename_p = sprintf('%s_%d.tle',Const_name,p);
web = 'https://celestrak.org/NORAD/elements';
url = sprintf('%s/gp.php?GROUP=%s&FORMAT=tle',web,Const_name);
cd save\tle\

% save TLE as 
filename_cc = sprintf('%s_cc.tle',Const_name);
if update_TLE == true % do it only once, or update every 2 hours!
    websave(filename_cc,url); % don't abuse -> server 403
end
try
    TLE_file = readlines(filename_cc) ; 
catch
    cd ../..
    error('Please, activate "update TLE online"')
end
    writelines(TLE_file,filename)

TLE = readlines(filename); clear TLE_file
N_tot = floor(length(TLE)/3); % total number of satellites available

% 2BP correction try: d(n), d2(n), B* = 0
for k = 1:N_tot
    line1_0 = char(TLE(2+(k-1)*3,:));
    
    line1_0(34:43) = ' .00000000'; % d(n)/d(t)
    line1_0(45:52) = ' 00000+0';   % d2(n)/d2(t)
    line1_0(54:61) = ' 00000+0';   % B*

    TLE(2+(k-1)*3,:) = string(line1_0);
end


if strcmpi(Max_sat,'All') % case insensitive string comparison
    Max_sat = N_tot;
end
N_pg  = ceil(N_tot/Max_sat);  % number of TLE pages

if p < N_pg
    writelines(TLE(1+(p-1)*3*Max_sat:p*3*Max_sat),filename_p)
elseif p == N_pg
    writelines(TLE(1+(p-1)*3*Max_sat:end),filename_p)
elseif p > N_pg
    error('try lower page number')
end


TLE_p = readlines(filename_p);
N = floor(length(TLE_p)/3); % Max_sat
muE = 398600.4418; % [km^3/s^2] 
RE = 6378.1366; % [km] Earth equatorial radius

for k = 1:N
    line0 = string(TLE_p(1+(k-1)*3,:));
    line1 = char(TLE_p(2+(k-1)*3,:));
    line2 = char(TLE_p(3+(k-1)*3,:));

    YY  = str2num(line1(19:20));
    DDD = str2num(line1(21:32));
    e   = str2num(['0.',line2(27:33)]); % eccentricity
    i   = str2num(line2(9:16));  % [deg] inclination 
    om  = str2num(line2(35:42)); % [deg] argument of periapsis
    OM  = str2num(line2(18:25)); % [deg] right ascension of the ascending node
    M   = str2num(line2(44:51)); % [deg] mean anomaly at epoch 
    n   = str2num(line2(53:63)); % [rev/day] mean motion
    T   = 1/n *24*60; % [min] orbital period
    h   = (muE*(T*60/(2*pi))^2)^(1/3) - RE; % [km] altitude (for circular orbits - calculated from semimajor-axis, assume R = RE_equator)
    if YY >= 57 % (first satellite launched in 1957)
        Epoch = datetime(sprintf('19%d-1-1',YY)) + DDD - 1;
    else
        Epoch = datetime(sprintf('20%d-1-1',YY)) + DDD - 1;
    end
    Name(k,:) = line0;
    DATA(p,1).Epoch(k,:) = Epoch;
    DATA(p,1).n(k,:) = n;
    DATA(p,1).T(k,:) = T;
    DATA(p,1).h(k,:) = h;
    DATA(p,1).e(k,:) = e;
    DATA(p,1).i(k,:) = i;
    DATA(p,1).om(k,:) = om;
    DATA(p,1).OM(k,:) = OM;
    DATA(p,1).M(k,:) = M;

    if h <= htr % assume satellites fall under this value
        DATA(p,1).off(k,:) = 1;
    else
        DATA(p,1).off(k,:) = 0;
    end
end
DATA(p,1).T_unique = unique(round(DATA(p,1).T,1));
DATA(p,1).h_unique = unique(round(DATA(p,1).h));
DATA(p,1).e_unique = unique(round(DATA(p,1).e,4));
DATA(p,1).i_unique = unique(round(DATA(p,1).i));
DATA(p,1).om_unique = unique(round(DATA(p,1).om));
DATA(p,1).OM_unique = unique(round(DATA(p,1).OM));
% assignin('base','DATA',DATA)

nominal = find(DATA(p,1).off == 0); % nominal satellites indices
off_sat = find(DATA(p,1).off == 1);
N = length(nominal);
N_off = length(off_sat);
if N_off > 0
    N_tot = length(DATA(p).off);
    index = repmat(1:N_tot,3,1);
    index = index(:);
    [INDEX,~] = find(index == nominal');
    writelines(TLE_p(INDEX),filename_p)
end

% Name = strrep(Name,' ',''); 
Name = deblank(Name); % removes (trailing) spaces
DATA(p,1).Name = Name;
DATA(p,1).num = ceil(DATA(p,1).T*60/SampleTime);
% mode / median / mean
NUM = mode(DATA.num); % most frequent value 
% T = mean(DATA.T);
% T = mode(DATA.T);
T = max(DATA.T);
app.DATA = DATA;

EPOCH   = sort(DATA(p,1).Epoch(nominal));
Epoch_0 = EPOCH(1);

StartTime = Epoch_0;
StopTime  = Epoch_0 + Sim_Time; 

if SAT_search == true && Walker_const == false
    % Name -> sat number
%     SAT = find(strcmpi(Name,Sat_name));
    SAT = find(contains(Name,Sat_name,'IgnoreCase',true));
    DATA_SAT.Epoch = DATA.Epoch(SAT);
    DATA_SAT.n = DATA.n(SAT);
    DATA_SAT.T = DATA.T(SAT);
    DATA_SAT.h = DATA.h(SAT);
    DATA_SAT.e = DATA.e(SAT);
    DATA_SAT(p,1).i = DATA.i(SAT);
    DATA_SAT(p,1).om = DATA.om(SAT);
    DATA_SAT(p,1).OM = DATA.OM(SAT);
    DATA_SAT.M = DATA.M(SAT);
    DATA_SAT.off = DATA.off(SAT);
    DATA_SAT.T_unique = unique(round(DATA_SAT.T,1));
    DATA_SAT.h_unique = unique(round(DATA_SAT.h));
    DATA_SAT.e_unique = unique(round(DATA_SAT.e,4));
    DATA_SAT(p,1).i_unique = unique(round(DATA_SAT(p,1).i));
    DATA_SAT(p,1).om_unique = unique(round(DATA_SAT(p,1).om));
    DATA_SAT(p,1).OM_unique = unique(round(DATA_SAT(p,1).OM));
    DATA_SAT.Name = Name(SAT,:);
    DATA_SAT.num = DATA.num(SAT);
%   mode / median / mean
    NUM = mode(DATA_SAT.num); % most frequent value  
%     T = mean(DATA_SAT.T);
%     T = mode(DATA_SAT.T);
    T = max(DATA_SAT.T);
    app.DATA_SAT = DATA_SAT;
    nominal = find(DATA_SAT.off == 0);
    app.SatellitesNumberTextArea.Value = join(string(SAT(nominal)),',',1);
    N = length(nominal);
    if ismember(1,DATA_SAT.off)        
        N_off = length(DATA_SAT.off(DATA_SAT.off==1));
    end
end


fprintf('%8.2f seconds is TLE Retrieval time.\n',toc)

%% SIMULATION - setup
tic
T_max = max(DATA(p,1).T);

tleFile = filename_p;

if SAT_search == true && Walker_const == false
    StartTime = sort(DATA_SAT.Epoch(nominal));
    StartTime = StartTime(1);
    StopTime  = StartTime + Sim_Time; 

    filename_SAT = sprintf('%s_SAT.tle',Sat_name);
    i = 3*SAT(nominal)-2;
%     i = 3*SAT-2;

    for j = 1:length(i)
        TLE_(3*j-2:3*j) = TLE(i(j):i(j)+2);
    end
    writelines(TLE_,filename_SAT)

    tleFile = filename_SAT;
    Name = DATA.Name(SAT,:);
    T_max = max(DATA_SAT.T);
    N = floor(length(TLE_)/3); 
end
if Walker_const == false
    app.NofSatsTextArea.Value = string(N);
    app.offSatsTextArea.Value = string(N_off);
end

Scenario = satelliteScenario(StartTime,StartTime,SampleTime);

constellation = satellite(Scenario,tleFile,'OrbitPropagator','two-body-keplerian','Name',Name(nominal)); % recommended, sgp4/sdp4 propagation for first state

cd ../..

fprintf('%8.2f seconds is Constellation Scenario time.\n',toc)

elseif Walker_const == true
    StartTime = datetime(StartDate);
    StopTime = StartTime + Sim_Time;
end
fprintf(' StartTime: %s\n StopTime:  %s\n',StartTime,StopTime)

%% STORAGE - Simulation starts
    tic
        
    t_utc = (StartTime:seconds(30):StopTime)';
    if Walker_const == false
        [r,v] = states(constellation,StartTime,'CoordinateFrame','inertial');
        r = (permute(r,[3,1,2])); % [m]  nsat x 3 x ntime
        v = (permute(v,[3,1,2])); % [m/s]
        fprintf('%8.2f seconds is Matlab simulation time.\n',toc)
    end
    
Steps_tot = ceil(seconds(Sim_Time)/SampleTime+1e-8);

%% SIMULINK 2BP/LLA CYCLE BUILD
Sim_Time = seconds(Sim_Time);
T_prog = ceil(Sim_Time/(10000*SampleTime))*SampleTime; % check progress sample time
t_UTC = datevec(t_utc);
muE = 398600.4418; % [km^3/s^2] 

if Walker_const == false
    ECI__ = cell(1);
else
    ECI__ = cell(length(a),1);
    TT = TT_P .* P; % total number of satellites (Walker orbit)
end
for jj = 1:length(ECI__)
    if Walker_const == true
        missionSatellites.SemiMajorAxis = a(jj)*ones(TT(jj),1)*1e3; % [m]
        missionSatellites.Eccentricity  = e(jj)*ones(TT(jj),1);
        missionSatellites.Inclination   = i(jj)*ones(TT(jj),1); % [deg]
        missionSatellites.ArgOfPeriapsis= om(jj)*ones(TT(jj),1); % [deg]
        S(jj) = TT(jj)/P(jj);
        missionSatellites.RAAN         = sort(repmat(0:360/P(jj):360-1e-8,1,S(jj)))'; % [deg]
        missionSatellites.TrueAnomaly  = repmat(0:360/S(jj):360-1e-8,1,P(jj))'; % [deg]
        df(jj) = F(jj)*360/TT(jj); % [deg] relative angular shift between adjacent orbital planes
        for j = 1:P(jj)-1
            missionSatellites.TrueAnomaly(j*S(jj)+1:(j+1)*S(jj)) = missionSatellites.TrueAnomaly(j*S(jj)+1:(j+1)*S(jj)) + df(jj)*j;
        end
        T(jj) = 2*pi*sqrt(a(jj)^3/muE)/60; % [min] orbital period
        N = TT(jj);

        missionSatellites.ConstellationDefinition = table(missionSatellites.SemiMajorAxis, ...
            missionSatellites.Eccentricity, missionSatellites.Inclination, ...
            missionSatellites.RAAN, missionSatellites.ArgOfPeriapsis, ...
            missionSatellites.TrueAnomaly, ...
            'VariableNames', ["a (m)", "e", "i (deg)", "Ω (deg)", "ω (deg)", "θ (deg)"]);
        GUI_DATA.missionSatellites(jj) = missionSatellites;
    end

    if AUTO_period_mod == true
    
    % Chobotov method, page 415 (Orbital Mechanics, Third Edition), here ratio = Q
    ome = 0.250684454; % [deg/min] rotation rate of Earth
    S = T(jj)*ome; % [deg] Westward shift on each pass

    ratio(jj) = 360/S;
    [N_(jj),D(jj)] = rat(ratio(jj)); % D total number of sidereal days until GT repeats itself, 
                                     % N_ number of sat revolutions until repetition, 
        if strcmpi(free_memory,'max')
        else
            n1(jj) = 0;
            free_memory = str2double(free_memory);
            D_max = floor(free_memory/(N*3*8*3600*23.934/SampleTime)); % max D allowed by memory, single precision, 
            % UPDATE: SIMULINK WORKS IN DOUBLE PRECISION, IT NEEDS MORE MEMORY!
            % THEREFORE USE 8 INSTEAD OF 4
            while D(jj) > D_max
                n1(jj) = n1(jj)+1;
                ratio(jj) = round(ratio(jj)*10^(10-n1(jj)))/10^(10-n1(jj));
                [N_(jj),D(jj)] = rat(ratio(jj));
            end
        end
        T_LLAcycle = hours(D(jj)*23.934);
    end
    if Periodcycle_mode == true && AUTO_period_mod == false
        T_LLAcycle = TTT;

    elseif Periodcycle_mode == false
        T_LLAcycle = hours(hours(seconds(Sim_Time)));
    end
    if T_LLAcycle > seconds(Sim_Time)
        T_LLAcycle = hours(hours(seconds(Sim_Time)));
    end


    tic
    if ECEF_mode == 0
        rf = 'ECI';
    else
        rf = 'ECEF';
    end
    if Walker_const == false
        model_name = 'ECI_propagation.slx';
        fprintf(' %s 2BP TLE SIMULINK simulation:\n',rf)
    elseif Walker_const == true
        model_name = 'Walker_Delta.slx';
        fprintf(' %s 2BP Walker-Delta SIMULINK simulation:\n',rf)
        r = [0,0,0];
        v = [0,0,0];
    end
    load_system(model_name)
    wkspace = sprintf('%c',model_name(1:end-4));
    mdlWks = get_param(wkspace,'ModelWorkspace');
    clear(mdlWks)
    simIn = Simulink.SimulationInput(wkspace);
    set_param(wkspace,'SimulationMode','normal')

    Ti = "0";
    Tf = num2str(str2double(Ti) + seconds(T_LLAcycle));
    simIn = setModelParameter(simIn,'StartTime',Ti,'StopTime',Tf);
    simIn = setVariable(simIn,'r',r,'Workspace',wkspace);
    simIn = setVariable(simIn,'v',v,'Workspace',wkspace);
    simIn = setVariable(simIn,'SampleTime',SampleTime,'Workspace',wkspace);
    simIn = setVariable(simIn,'Sim_Time',Sim_Time,'Workspace',wkspace);
    simIn = setVariable(simIn,'T_prog',T_prog,'Workspace',wkspace);
    simIn = setVariable(simIn,'StartTime',StartTime,'Workspace',wkspace);
    if  Walker_const == true
        simIn = setVariable(simIn,'StartDate',StartDate,'Workspace',wkspace);
        simIn = setVariable(simIn,'missionSatellites',missionSatellites,'Workspace',wkspace);
    end

    applyToModel(simIn)

    if ECEF_mode == 0
        set_param([wkspace,'/Orbit Propagator Numerical (high precision)'],'Propagator','kepler','outportFrame','ICRF')
    else
        set_param([wkspace,'/Orbit Propagator Numerical (high precision)'],'Propagator','numerical','outportFrame','Fixed-frame')
    end

    % warning off aerospace:orbit:KeplerFailedToConverge
    out = sim(simIn);
    fprintf('\b\b\b\b\b\b\b\b%7.3f%%',100)
    fprintf('\n%8.2f seconds is SIMULINK %d time.\n',toc,jj)
    save_system(model_name)
    close_system(model_name)
    ECI__{jj,1} = out.ECI(:,1:3,:);
    if size(ECI__{jj,1},3) == 1
        ECI__{jj,1} = permute(ECI__{jj,1},[3 2 1]);
    end
end
if Walker_const == true
    T = seconds(minutes(T)); 
    app.T = T;
end
% CYCLING LLA
if ECEF_mode == 0
    t_utc0 = datevec(StartTime);

    LLA__ = cell(size(ECI__));
    for jj = 1:length(ECI__)
        fprintf(' ECI2LLA conversion:        ')
        tic
        n00 = 0;
        LLA_ = single(zeros(size(ECI__{jj,1})));
        sz3 = size(ECI__{jj,1},3);
        for k = 1:sz3
            LLA_(:,:,k) = ECI2LLA(ECI__{jj,1}(:,:,k)*1e-3,t_utc0,(k-1)*SampleTime); % [rad,rad,km]
            if abs(n00 - k) > T_prog % sampling steps for progression stamp, it would break if checked at every step
                n00 = k;
                fprintf('\b\b\b\b\b\b\b%6.2f%%',k/sz3*100)
            end
        end
        LLA__{jj,1} = LLA_;
        fprintf('\b\b\b\b\b\b\b%6.2f%%',100)
        fprintf('\n%8.2f seconds is LLA %d array time.\n',toc,jj)
    end
    % end

    sz = size(LLA__,1);
    if Periodcycle_mode == true
        tic
        LLA_ = [];
        for jj = 1:sz
            LLA_ = cat(1,LLA_,LLA__{jj});
        end
        fprintf('%8.2f seconds is LLA given Period time.\n',toc)
    end
end

%% MAIN LOOP

% fprintf(' Main loop:        ')
%     n00 = 0;
%     tic
%     for k = 1:Steps_tot
%         LLA = [];
%         for jj = 1:sz
%             sz3 = size(LLA__{jj,1},3);
%             index = mod(k,sz3);
%             I = find(~index);
%             index(I) = sz3(I);
%             LLA = cat(1,LLA,LLA__{jj}(:,:,index));
%         end
% 
% %         end
%         if abs(n00 - k) > T_prog % sampling steps for progression stamp, it would break if checked at every step
%             n00 = k;
%             fprintf('\b\b\b\b\b\b\b%6.2f%%',k/Steps_tot*100)
%         end
%     end
%     fprintf('\b\b\b\b\b\b\b%6.2f%%',100)
% fprintf('\n%8.2f seconds is LLA cycling time.\n',toc)





%% DATA to Base Workspace and GUI

app.ECI = ECI__;
if ECEF_mode == 0
    app.LLA = LLA__;
end

GUI_DATA.app = app;
if Walker_const == false
    GUI_DATA.Const_name = Const_name;
    GUI_DATA.Max_sat = Max_sat;
    GUI_DATA.p = p;
    GUI_DATA.StartDate = datevec(StartTime);
    GUI_DATA.Sim_Time = Sim_Time ;
    GUI_DATA.SampleTime = SampleTime;
    GUI_DATA.update_TLE = update_TLE;
    GUI_DATA.SAT_search = SAT_search;
    GUI_DATA.Sat_name = Sat_name;
    GUI_DATA.DATA = DATA;
    DATA.Name = DATA.Name(nominal);
    if SAT_search == true
        GUI_DATA.DATA_SAT = DATA_SAT;
        DATA_SAT.Name = DATA_SAT.Name(nominal);
    end
else
    GUI_DATA.Const_name = 'Walker-Delta';
    GUI_DATA.i = i;
    GUI_DATA.TT = TT;
    GUI_DATA.P = P;
    GUI_DATA.F = F;
    GUI_DATA.a = a;
    GUI_DATA.e = e;
    GUI_DATA.om = om;
    GUI_DATA.StartDate = StartDate;
    GUI_DATA.Periodcycle_mode = Periodcycle_mode;
    GUI_DATA.TTT = TTT;
end

GUI_DATA.savevarname = savevarname;
GUI_DATA.Sat_to_plot = Sat_to_plot;
GUI_DATA.free_memory = free_memory;         
GUI_DATA.folderpath = folderpath;
GUI_DATA.GS = GS;

assignin('base','GUI_DATA',GUI_DATA)
assignin('base','coord_es',coord_es)
cd save\
if ECEF_mode == 0
    assignin('base',"ECI",ECI__)
    if Periodcycle_mode == false
        assignin('base',savevarname,LLA__)
        save([savevarname,'.mat'],'LLA__','-mat','-v7.3')
    else
        assignin('base',savevarname,LLA_)
        save([savevarname,'.mat'],'LLA_','-mat','-v7.3')
    end
else
    assignin('base',"ECEF",ECI__)
end
cd ../

t = toc(tic0);
fprintf('%8.2f seconds is Total time.\n',t)

app.ExectimeTextArea.Value = sprintf('%.1f seconds',t);


%% PLOTS
close all

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


if length(ECI__) == 1 %|| Walker_const == false
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
    grid(axes,'on')
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
exportgraphics(gca,sprintf('%s.png',frm))
cd ../..

% close figures, plots have been saved and are available in app only
close(1)
close(2)

%% FUNCTIONS

% function y = ode_2bp(~,x)
% 
% x = reshape(x,6,[])';
% 
% muE = 3.986004354360959e+05; % [km^3/s^2] from spice  
% r = x(:,1:3);
% v = x(:,4:6);
% % y = zeros(6,1);
% y(:,1:3) = v;
% y(:,4:6) = -muE*r./(vecnorm(r,2,2).^3); % 2bp dynamics
% 
% y = reshape(y',[],1);
% 
% end
% 
% %--------------------------------------------------------------------------
% 
% function [r,v] = orbit_2bp_propagation(x0,t)
% % x0 is n x 6, n satellites, and has to be provided to ode113 
% % as a vector
% 
% x0 = 1e-3*reshape(x0',[],1); % [m -> km]
% 
% options = odeset('RelTol',1e-13,'AbsTol',1e-14);
% [~,x] = ode113(@ode_2bp,t,x0,options);
% 
% x = 1e3*reshape(x(end,:),6,[])'; % [km -> m]
% r = x(:,1:3);  % rows are different satellites, single time
% v = x(:,4:6);
% 
% end


%--------------------------------------------------------------------------

function LLA = ECI2LLA(ECI,t_utc0,t)
% INPUTs:
%   ECI    : [km] x,y,z, nx3 (n satellites)
% OUTPUTs:  mode I
%       t_utc0 : 1x6 initial time vector [YYYY,MM,DD,hh,mm,ss], 1901 <= YYYY <= 2099
%       t      : [s] vector time step (used as scalar)
%           mode II
%       t_utc0 : scalar, tG0 Greenwich initial time Longitude [rad] -> used
%       t      : [s] vector time step (used as scalar)


% OUTPUTs:
%   LLA    : geodetic Lat,Lon,Alt [rad,rad,km]

if nargin == 2
    t = 0;
end
if size(t_utc0,2) == 6
    % Curtis (4th Edition) method for Greenwich Sidereal Time pp. 249-252,
    % checked with ECEF from Simulink Orbit Propagator block
    J0 = 367*t_utc0(1)-fix(7/4*(t_utc0(1)+fix((t_utc0(2)+9)/12)))+fix(275/9*t_utc0(2))+t_utc0(3)+1721013.5;
    T0 = (J0-2451545)/36525;
    lonG0 = mod(100.4606184 + 36000.77004 * T0 + 0.000387933 * T0^2 - 2.583e-8 * T0^3,360); % [deg]
    lonG = mod(lonG0 + 360.98564724 * (t_utc0(4) + t_utc0(5)/60 + t_utc0(6)/3600)/24,360); % [deg]
    tG0 = lonG * pi/180; % [rad]
    % tG0 = 0;
   
elseif size(t_utc0,2) == 1  
    tG0 = t_utc0; 
end

omE = pi/180/3600*15.04; % Earth radial velocity [rad/s]

lonG = zeros(size(ECI,1));
alpha = zeros(size(ECI,1),length(t));
lon = alpha; lat = alpha; lat_gd = lat; h = lat;
% delta = lonG;

    lonG = mod(tG0 + t*omE, 2*pi);  % [rad]
    x = ECI(:,1); 
    y = ECI(:,2);
    z = ECI(:,3);
    rn = vecnorm(ECI,2,2);   % [km]

    lat(:) = asin(z./rn);  % [rad]  Declination

    alpha(:) = atan2(y,x);  % [rad]  Right Ascension

    lon(:) = alpha(:) - lonG;

 I = find(lon(:) > pi);
 lon(I)  = lon(I) - 2*pi;

 I = find(lon(:) < -pi);
 lon(I)  = lon(I) + 2*pi;

% GEODETIC LATITUDE - LONGITUDE - ALTITUDE
% geodetic Longitude = geocentric Longitude
%     [lat_gd(:,k),h(:,k)] = geoc2geod(lat(:,k)*180/pi,rn*1e3);

%     [lat_gd(:),h(:)] = geoc2geod(lat(:)*180/pi,rn*1e3,'WGS84');
%     lat_gd(:) = lat_gd(:)*pi/180; % [rad]
%     h(:) = h(:)*1e-3; % [km]  Ellipsoidal / Geodetic Altitude
 
% end
h = rn - 6371; % geometric altitude over mean Earth radius
LLA = single([lat,lon,h]);

% LLA = [lat_gd,lon,h];
% LLA = single([lat_gd,lon,h]);

end

%--------------------------------------------------------------------------

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

xticks(axes,-180:30:180)
yticks(axes,-90:30:90)
set(axes,'YDir','normal')
hold(axes,'on')
grid(axes,'on')
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