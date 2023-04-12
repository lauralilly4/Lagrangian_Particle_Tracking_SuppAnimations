%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Code for Supplemental Animations to:
%%%     Lilly et al. Using a Lagrangian particle tracking model to evaluate 
%%%     impacts of El Niño-related advection on euphausiids in the southern 
%%%     California Current System
%%% Author: Laura E. Lilly
%%% Updated: 27 Jan 2022

% CODING STEP 3: Run particle backtracks


%% ======= Step 0: Parameters and files =======
yrin = cutyr;
sppin = spptype;
mdst = [3,31];
mded = [12,1];
tStep = 0.01;
zdp = 1:15;
conc1 = 100;

% Load files
load FMT.mat;
load MIT_CCS_LATLON_500m.mat;
load grid.mat;
[nxc nyc nzc] = size(hFacC);
nt = 3743;
time = rdslice(['TIME_CCS_2007_2020.bin'],[nt 1],1, fmt, Ieee);

% %%% LAT/LON BOUNDARIES %%%
[lon0,lon0loc] = min(XC(:)); % Degrees E, btw
[lon1,lon1loc] = max(XC(:));
[lat0,lat0loc] = min(YC(:));
[lat1,lat1loc] = max(YC(:));

% Set lon & lat limits to be full region
xrmin = lon0;
xrmax = lon1;
yrmin = lat0;
yrmax = lat1;
plons1 = (xrmin:(xrmax-xrmin)/conc1:xrmax);
plats1 = (yrmin:(yrmax-yrmin)/conc1:yrmax);
[xinit,yinit] = meshgrid(plons1,plats1);
zinit = xinit * 0;  
    
    
%% ======= Time matchup to get model chunk =======
dstart = datenum(datetime(yrin,mdst(1),mdst(2),12,0,0));

if mded(1) == 1
    dend = datenum(datetime(yrin,mded(1),mded(2),12,0,0));
elseif mded(1) == 12
    dend = datenum(datetime(yrin-1,mded(1),mded(2),12,0,0));
end

tsid = find(time == dstart);
teid = find(time == dend);
tchunk = time(tsid:-1:teid); % Vector of DATES of time-chunk - REVERSE if indices inputted are from Mar 31 -> Jan 1
tidx = tsid:-1:teid; % Vector of INDICES of time-chunk
tneg = (-1)*(tchunk); % Time-chunk is already in reverse order because inputted 
                      % indices were that way -> now just multiply by -1 to make numbers
                      % in ascending order (but all neg)

                      
%% ======= Step 1: Particle Movements Loop -> Threshold #1 (>80%) =======
% Allocate memory for particle time-tracks
xp = NaN(length(xinit),length(xinit),length(tchunk)); % x, y, timestep
yp = NaN(length(xinit),length(xinit),length(tchunk));
zp = NaN(length(xinit),length(xinit),length(tchunk));

% Allocate memory for 2 adjacent u/v slices
uslc = NaN(nxc,nyc,2);
vslc = NaN(nxc,nyc,2);
wslc = NaN(nxc,nyc,2);

inc = tsid+1; % start with 'idx1+1' so I can immediately do +1

for ii = 1:(length(tchunk)-1)

    inc = inc-1;

    if(inc == tsid)

        xp(:,:,ii) = xinit; 
        yp(:,:,ii) = yinit;
        zp(:,:,ii) = zinit;

        % Get first and second U&V slices 
        ucut1 = rdslice(['U_MIT_CCS_2007_2020_fulldomain.bin'],[nxc nyc nzc],inc, fmt, Ieee);
        vcut1 = rdslice(['V_MIT_CCS_2007_2020_fulldomain.bin'],[nxc nyc nzc],inc, fmt, Ieee);
        uslc(:,:,1) = (-1)*nanmean(ucut1(:,:,zdp),3);
        vslc(:,:,1) = (-1)*nanmean(vcut1(:,:,zdp),3);
        wslc(:,:,1) = 0 * uslc(:,:,1);

        ucut2 = rdslice(['U_MIT_CCS_2007_2020_fulldomain.bin'],[nxc nyc nzc],inc-1, fmt, Ieee);
        vcut2 = rdslice(['V_MIT_CCS_2007_2020_fulldomain.bin'],[nxc nyc nzc],inc-1, fmt, Ieee);
        uslc(:,:,2) = (-1)*nanmean(ucut2(:,:,zdp),3);
        vslc(:,:,2) = (-1)*nanmean(vcut2(:,:,zdp),3);
        wslc(:,:,2) = 0 * uslc(:,:,2);    

    else
        % Don't need the setup particle position stuff
        ucut1 = rdslice(['U_MIT_CCS_2007_2020_fulldomain.bin'],[nxc nyc nzc],inc, fmt, Ieee);
        vcut1 = rdslice(['V_MIT_CCS_2007_2020_fulldomain.bin'],[nxc nyc nzc],inc, fmt, Ieee);
        uslc(:,:,1) = (-1)*nanmean(ucut1(:,:,zdp),3);
        vslc(:,:,1) = (-1)*nanmean(vcut1(:,:,zdp),3);
        wslc(:,:,1) = 0 * uslc(:,:,1);

        ucut2 = rdslice(['U_MIT_CCS_2007_2020_fulldomain.bin'],[nxc nyc nzc],inc-1, fmt, Ieee);
        vcut2 = rdslice(['V_MIT_CCS_2007_2020_fulldomain.bin'],[nxc nyc nzc],inc-1, fmt, Ieee);
        uslc(:,:,2) = (-1)*nanmean(ucut2(:,:,zdp),3);
        vslc(:,:,2) = (-1)*nanmean(vcut2(:,:,zdp),3);
        wslc(:,:,2) = 0 * uslc(:,:,2);    
    end
   
    % We now have 2 time slices of data (current time, time+1)
    % Integrate using an RK4 scheme and interpolation


    %% 2. Compute gridded interpolants for pairs
    x = XC;
    y = YC;
    z = ZC;
    ti = tneg(ii:ii+1);

    [gi.x,gi.y,gi.ti] = ndgrid(x(:,1),y(1,:),ti);
    gi.F_u = griddedInterpolant(gi.x,gi.y,gi.ti,uslc);
    gi.F_v = griddedInterpolant(gi.x,gi.y,gi.ti,vslc);
    gi.F_w = griddedInterpolant(gi.x,gi.y,gi.ti,wslc);


    %% 3. Compute particle tracks
    %%%% tStep is h in flow_map_3d - therefore, you should make it as coarse as
    %%%% possible without affecting your tracks. The interpolation scheme
    %%%% works from the output, but refines the space & time grid to better
    %%%% advect particles.
    tStart = tneg(ii); % Initial time ***DIFFERENT than 'PrtTrcks'
    tFinal = tneg(ii+1); % Final time

    % Use previous scheme (Michael Allshouse): 2D interpolation
    [xTrajectory, yTrajectory, zTrajectory] = ...
        ch3_flow2d(xp(:,:,ii),yp(:,:,ii),zp(:,:,ii),...
        tStart,tFinal,tStep,gi); % the integration

    % Calculate the displacement of each particle over this step
    dxp = xTrajectory(:,:,end) - xTrajectory(:,:,1); 
    dyp = yTrajectory(:,:,end) - yTrajectory(:,:,1);
    dzp = zTrajectory(:,:,end) - zTrajectory(:,:,1);

    % Calculate the new particle position
    xp(:,:,ii+1) = xp(:,:,ii) + dxp;
    yp(:,:,ii+1) = yp(:,:,ii) + dyp;
    zp(:,:,ii+1) = zp(:,:,ii) + dzp; 

end 


% Total particle displacement
dxp = xp(:,:,end)-xp(:,:,1);
dyp = yp(:,:,end)-yp(:,:,1);
dzp = zp(:,:,end)-zp(:,:,1);

% Vectorize xp and yp displacements
xpt = NaN(length(xp(:,1,1)),length(xp(1,1,:)));
ypt = NaN(length(xp(:,1,1)),length(xp(1,1,:)));

for xt = 1:length(xp(1,1,:))
    xd = diag(xp(:,:,xt));
    yd = diag(yp(:,:,xt));
    xpt(:,xt) = xd;
    ypt(:,xt) = yd;
end


% Squeeze 3D particle positions -> 2D
xpqz1 = NaN(length(xp(:,1,1))^2,length(xp(1,1,:)));
ypqz1 = NaN(length(yp(:,1,1))^2,length(yp(1,1,:)));

for xq = 1:length(xp(1,1,:))
    sqx = reshape(xp(:,:,xq),[1,length(xp(:,1,1))^2]);
    sqy = reshape(yp(:,:,xq),[1,length(yp(:,1,1))^2]);
    xpqz1(:,xq) = sqx;
    ypqz1(:,xq) = sqy;
end

% Get indices of particle locations within polygon
thrshidx1 = find(inpolygon(xpqz1(:,1),ypqz1(:,1),lonsmin1(thrsh1bnd),latsmin1(thrsh1bnd)));


% Get indices of particle locations in polygon -> within >50% but NOT in
% >80% (inner parts from 'thrshidx1')
thrshidx2all = find(inpolygon(xpqz1(:,1),ypqz1(:,1),lonsmin2(thrsh2bnd),latsmin2(thrsh2bnd)));
xpqz2sub = xpqz1(thrshidx2all,:);
ypqz2sub = ypqz1(thrshidx2all,:);
thrshidx2 = find(~inpolygon(xpqz2sub(:,1),ypqz2sub(:,1),lonsmin1(thrsh1bnd),latsmin1(thrsh1bnd)));



%% ======= PLOT: animation of backtracks - SUBSET region =======
% Save as .avi (not .gif)
clear F
clear frame
clear newVid

szstep = 4;
h = figure(3)

for pq = 1:length(xpqz1(1,:))

    % Plot >80% region first
    scatter(xpqz1(thrshidx1,pq),ypqz1(thrshidx1,pq),szstep+8,'o','b','filled'); % Thresh #1 - end

    hold on
    % Then plot 50-80% region
    scatter(xpqz2sub(thrshidx2,pq),ypqz2sub(thrshidx2,pq),szstep+6,'d','filled','MarkerFaceColor',[0 0.7 0.7],'MarkerEdgeColor',[0 0.7 0.7]);

    plot(coastlonall+360,coastlatall,'k','LineWidth',2);
    hold off
    
    title([sppin,' - ',datestr(time(tidx(pq)),'mm-dd-yyyy')],'FontSize',24);
    xlim([min(XC(:)),max(XC(:))]);
    ylim([min(YC(:)),max(YC(:))]);
    set(gca,'dataaspectratio',[1 cos(pi/180*(min(YC(:))+max(YC(:)))/2) 1]);
    
    F(pq) = getframe(gcf);
end



% Create a video object to fill
newVid = VideoWriter(strcat(spptype,'_Mar',datestr(time(tidx(1)),'yyyy'),'_Dec',datestr(time(tidx(1))-365,'yyyy'),'_hindcast.avi'));
newVid.FrameRate = 2;
newVid.Quality = 100;
open(newVid)

% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(newVid, frame.cdata);
end
% close the writer object <- VERY IMPORTANT!!!!!
close(newVid);
