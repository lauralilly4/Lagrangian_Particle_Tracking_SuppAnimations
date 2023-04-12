%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Code for Supplemental Animations to:
%%%     Lilly et al. Using a Lagrangian particle tracking model to evaluate 
%%%     impacts of El Niño-related advection on euphausiids in the southern 
%%%     California Current System
%%% Author: Laura E. Lilly
%%% Updated: 27 Jan 2022

% CODING STEP 2: Load euphausiid species, objectively map distribution,
% calculate 50% and 80% threshold boundaries


%% ======= Step 1: File import =======
spptype = input('Species?     '); % 
zoopfl = importdata(strcat(spptype,'_CSB_nt_tot_unt_bm.csv'));


%% ======= Step 2: Assign variables and get *ONE* value/station/spring =======
% % Lat, lon, abundance, biomass
lat = zoopfl.data(:,1);
lon = zoopfl.data(:,2);
bm = zoopfl.data(:,3);

% Dates, Line/Stn info
zoopyrs = [];
zoopmons = [];
zoopdys = [];
lines = [];
stns = [];
latlonsbms = [];


for y=2:length(zoopfl.textdata)
    
    % Check for and skip stations with extra letters ('E','U','C')
    if length(zoopfl.textdata{y,7}) > 2 &&  zoopfl.textdata{y,7}(3) == 'E'
       continue
    elseif zoopfl.textdata{y,7}(1) == 'U'
       continue
    elseif zoopfl.textdata{y,7}(1) == 'C'
       continue
    else
        yr = str2num(zoopfl.textdata{y,5}); % Years: convert 'cell' to number
        mon = str2num(zoopfl.textdata{y,5}); % Years: convert 'cell' to number  
        dy = str2num(zoopfl.textdata{y,4});
        lin = str2num(zoopfl.textdata{y,6});
        stn = str2num(zoopfl.textdata{y,7});
        latlonbm = [zoopfl.data(y-1,:)];

        zoopyrs = [zoopyrs; yr];
        zoopmons = [zoopmons; mon];
        zoopdys = [zoopdys; dy];
        lines = [lines; lin];
        stns = [stns; stn];
        latlonsbms = [latlonsbms;latlonbm];
    end
end     


% Array #1: UNSORTED Year, Month, Day, Line, Stn, Lat, Lon, Biomass
zoopdatesall = [zoopyrs,zoopmons,zoopdys,lines,stns,latlonsbms];
% Array #2 - SORT by Year, Mon, Line, Stn
zoopsort = sortrows(sortrows(sortrows(sortrows(zoopdatesall,5),4),2),1);

% Convert untransformed data -> log-transformed
zdatarl = zoopsort(~isnan(zoopsort(:,8)),:); % First, get only 'real' (non-NaN) values
znumlg = log10(zdatarl(:,8)+1); % Convert biomasses -> log-transform
zdatalg = [zdatarl(:,1:7),znumlg]; % Recombine log-transformed data w/ other variables

zdatarl = zoopsort(~isnan(zoopsort(:,8)),:); % First, get only 'real' (non-NaN) values
znumlg = log10(zdatarl(:,8)+1); % Convert biomasses -> log-transform
zdatalg = [zdatarl(:,1:7),znumlg]; % Recombine log-transformed data w/ other variables


%% One value per station per spring (in case multiple cruises and/or samples)
dcutin = zdatalg;
zcol = 8;

yrsunq = unique(dcutin(:,1)); % Get list of all years for for-loop
allzoopsubs = []; % Create empty array for each year's subset data (have to stack vertically because they will have different nos. rows
allcutarrs = [];
allmeanarrs = [];
allzoopmeans = [];

for yr=1:length(yrsunq)
    yridx = find(dcutin(:,1) == yrsunq(yr));
    yrsub = dcutin(yridx,:);

    mocutunq = unique(yrsub(:,2));

    % Get only ONE spring month of data (plus additional *UNIQUE* Ln/Stn combos from other months): 
    if length(mocutunq) == 1 % If only one month represented, just save all those rows
        zmocut4 = yrsub; 

    elseif length(mocutunq) > 1
        if sum(ismember(mocutunq,4)) > 0 % Check first if there are samples from Month = 4 (arbitrary priority)
            zmocut1 = yrsub(yrsub(:,2)==4,:);

            % Check for unique Ln/Stn samples in Month = 3
            if sum(ismember(mocutunq,3))> 0
                zmo3 = yrsub(yrsub(:,2)==3,:);
                zmo3(ismember([zmo3(:,4),zmo3(:,5)],[zmocut1(:,4),zmocut1(:,5)],'rows'),:) = [];
                zmocut2 = [zmocut1;zmo3];

                % Check for unique Ln/Stn samples in Month = 5 *in addition* to
                % Months 4,3
                if sum(ismember(mocutunq,5))> 0
                    zmo5 = yrsub(yrsub(:,2)==5,:);
                    zmo5(ismember([zmo5(:,4),zmo5(:,5)],[zmocut2(:,4),zmocut2(:,5)],'rows'),:) = [];
                    zmocut3 = [zmocut2;zmo5];

                    % Check for unique Ln/Stn samples in Month = 2 *in addition* to
                    % Months 4,3,5
                    if sum(ismember(mocutunq,2))> 0
                        zmo2 = yrsub(yrsub(:,2)==2,:);
                        zmo2(ismember([zmo2(:,4),zmo2(:,5)],[zmocut3(:,4),zmocut3(:,5)],'rows'),:) = [];
                        zmocut4 = [zmocut3;zmo2];
                    else
                        zmocut4 = zmocut3;
                    end

                elseif sum(ismember(mocutunq,2))> 0
                    zmo2 = yrsub(yrsub(:,2)==2,:);
                    zmo2(ismember([zmo2(:,4),zmo2(:,5)],[zmocut2(:,4),zmocut2(:,5)],'rows'),:) = [];
                    zmocut4 = [zmocut2;zmo2];
                else
                    zmocut4 = zmocut2;
                end

            elseif sum(ismember(mocutunq,5))> 0
                    zmo5 = yrsub(yrsub(:,2)==5,:);
                    zmo5(ismember([zmo5(:,4),zmo5(:,5)],[zmocut1(:,4),zmocut1(:,5)],'rows'),:) = [];
                    zmocut3 = [zmocut1;zmo5];

                    if sum(ismember(mocutunq,2))> 0
                        zmo2 = yrsub(yrsub(:,2)==2,:);
                        zmo2(ismember([zmo2(:,4),zmo2(:,5)],[zmocut3(:,4),zmocut3(:,5)],'rows'),:) = [];
                        zmocut4 = [zmocut3;zmo2];
                    else
                        zmocut4 = zmocut3;
                    end

            elseif sum(ismember(mocutunq,2))> 0
                    zmo2 = yrsub(yrsub(:,2)==5,:);
                    zmo2(ismember([zmo2(:,4),zmo2(:,5)],[zmocut1(:,4),zmocut1(:,5)],'rows'),:) = [];
                    zmocut4 = [zmocut1;zmo2];       
            else
                zmocut4 = zmocut1;
            end

        elseif sum(ismember(mocutunq,3))> 0
            zmocut1 = yrsub(yrsub(:,2)==3,:);

            if sum(ismember(mocutunq,5))> 0
                zmo5 = yrsub(yrsub(:,2)==5,:);
                zmo5(ismember([zmo5(:,4),zmo5(:,5)],[zmocut1(:,4),zmocut1(:,5)],'rows'),:) = [];
                zmocut2 = [zmocut1;zmo5];

                % Check for unique Ln/Stn samples in Month = 2 *in addition* to
                % Months 4,3,5
                if sum(ismember(mocutunq,2))> 0
                    zmo2 = yrsub(yrsub(:,2)==2,:);
                    zmo2(ismember([zmo2(:,4),zmo2(:,5)],[zmocut2(:,4),zmocut2(:,5)],'rows'),:) = [];
                    zmocut4 = [zmocut2;zmo2];
                else
                    zmocut4 = zmocut2;
                end

            elseif sum(ismember(mocutunq,2))> 0
                zmo2 = yrsub(yrsub(:,2)==5,:);
                zmo2(ismember([zmo2(:,4),zmo2(:,5)],[zmocut1(:,4),zmocut1(:,5)],'rows'),:) = [];
                zmocut4 = [zmocut1;zmo2];
            else
                zmocut4 = zmocut1;
            end

        elseif sum(ismember(mocutunq,5))> 0
            zmocut1 = yrsub(yrsub(:,2)==5,:);

            if sum(ismember(mocutunq,2))> 0
                zmo2 = yrsub(yrsub(:,2)==2,:);
                zmo2(ismember([zmo2(:,4),zmo2(:,5)],[zmocut1(:,4),zmocut1(:,5)],'rows'),:) = [];
                zmocut4 = [zmocut1;zmo2];
            else
                zmocut4 = zmocut1;
            end
        end
    end

    % If multiple entries exist for a ln/stn combo *from one month*, take only
    % Entry #1 (arbitrary -> gotta do it)
    [~,unqidx] = unique([zmocut4(:,4),zmocut4(:,5)],'rows');
    
    allstnavgs = [];
    for u=1:length(unqidx)
        allstns = zmocut4(zmocut4(:,4) == zmocut4(unqidx(u),4) & zmocut4(:,5) == zmocut4(unqidx(u),5),:);
        if length(allstns) > 1
            if length(allstns(1,:)) == 9
                stnavg = [allstns(1,1:7),nanmean(allstns(:,zcol)),nanmean(allstns(:,9))];
            elseif length(allstns(1,:)) == 8
                stnavg = [allstns(1,1:7),nanmean(allstns(:,zcol))];
            else
                continue
            end
        else
            stnavg = allstns;
        end
        allstnavgs = [allstnavgs;stnavg];
    end
     
    zoopcutunq = zmocut4(unqidx,:);  
    lnstncut = zoopcutunq(:,4)+zoopcutunq(:,5)*0.0001; % Combine line & stn into one value
    cutarr = horzcat(yrsunq(yr)*ones(length(zoopcutunq(:,1)),1),lnstncut,zoopcutunq(:,zcol:end));
    cutmean = horzcat(yrsunq(yr)*ones(length(zoopcutunq(:,1)),1),lnstncut,allstnavgs(:,end));
    allzoopsubs = [allzoopsubs; zoopcutunq]; % Save all yearly cuts -> ONE (arbitrary) from multiple samples
    allcutarrs = [allcutarrs; cutarr]; % Save each year's 'array' for Pivot Table
    allmeanarrs = [allmeanarrs; cutmean];
    
    % allzoopmeans: Year, Mon, Day, Line, Stn, Lat, Lon, Abd, Pred_abd
    allzoopmeans = [allzoopmeans;allstnavgs]; % All yearly cuts -> AVGed for multiple samples
end


%% ======== Step 3: Subset out single year cut =======
zdatavec = allzoopmeans;
zdata = zdatavec(:,zcol);
zlats = zdatavec(:,6);
zlons = zdatavec(:,7);

cutyr = input('Input year?  ');
entp = num2str(cutyr);
zdatasub = zdatavec(zdatavec(:,1)==cutyr,:); % CUT OUT: only data for that year

zlatcut = zdatasub(:,6);
zloncut = zdatasub(:,7); 
zdcut = zdatasub(:,zcol);


%% ======= Step 4: Calculate objective map parameters =======
xcors = 160;
ycors = 190; 
phi = -60; 
err = 0.1; 

xcor2s = xcors*xcors; 
ycor2s = ycors*ycors;
cosphi = cos(phi/180*pi);
sinphi = sin(phi/180*pi);

d = zdcut; % Input zoop dataset
zlatsub = zlatcut;
zlonsub = zloncut;

nxs = 100; % Number of grid points (x and y directions)
nys = 101;
zlons(zlons > 0) = -(zlons(zlons > 0)); % Change erroneous positive longitudes to negative
% Put in actual limits from full CalCOFI grid
alon1s = -126.717; 
alon2s = -112.117;
alat1s = 24.15;
alat2s = 37.967;

tokm = 111.324; % Conversion distance: lat -> km 
tokmt = tokm*cos(pi/180*(alat1s+alat2s)/2); % Conversion distance: lon -> km
xls = 0.0;
yls = 0.0; 
xrs = tokmt*(alon2s-alon1s); % Conversion of total x distance to km
yrs = tokm*(alat2s-alat1s); % Conversion of total y distance to km
xgs = linspace(xls,xrs,nxs);
ygs = linspace(yls,yrs,nys);
nfunc = 3;

% Convert distances to km and put data in vector d
xs = tokmt*(zlonsub-alon1s);
ys = tokm*(zlatsub-alat1s);
nfact = length(d);
disp(sprintf('%5i rows',nfact))

% %% ==== Transformation matrices ====
% Part 1: Calculate the data-data covariance matrix Eplus
[X1S,X2S] = meshgrid(xs);
[Y1S,Y2S] = meshgrid(ys);
dxs = X1S-X2S;
dys = Y1S-Y2S;
dxrs = dxs*cosphi-dys*sinphi;
dyrs = dxs*sinphi+dys*cosphi;
E = exp(-dxrs.^2/xcor2s-dyrs.^2/ycor2s);
Eplus = E+err*eye(nfact);

% Part 2: Calculate the data-grid points covariance matrix C.
C = zeros(nfact,nxs*nys);
[YS,XS] = meshgrid(ygs,xgs);
for n = 1:nfact
   dx2s = XS(:)'-xs(n);
   dy2s = YS(:)'-ys(n);
   dxr2s = dx2s*cosphi-dy2s*sinphi;
   dyr2s = dx2s*sinphi+dy2s*cosphi;
   C(n,:) = exp(-dxr2s.^2/xcor2s-dyr2s.^2/ycor2s);
end

% Part 3: Calculate the matrix A to turn data into a grid
m = 1:nfunc;
F = fxy(m,xs,ys);
EF = Eplus\F;
W = EF/(F'*EF);
EC = Eplus\C;
f = fxy(m,XS(:),YS(:))';
A = W*f+(eye(nfact)-W*F')*EC;

% Part 4: Calculate grids
% Data grid
grid = zeros(nxs,nys);
grid(:) = A'*d;

% Error grid
ergrid = zeros(nxs,nys);
EC = Eplus*A;
for n = 1:nxs*nys
   ergrid(n) = A(:,n)'*(EC(:,n)-2*C(:,n))+1;
end



%% ======= Step 5: Mask land & negatives -> PLOT Obj Map =======
thresh = 0.31;
gin = grid;

% Change negative values -> ZERO
nangrid = zeros(size(gin));
for g=1:length(gin(:,1))
    for r=1:length(gin(1,:))
        if gin(g,r) < 0;
            gridcut = 0;
        else
            gridcut = gin(g,r);
        end
        nangrid(g,r) = gridcut;
    end
end

pltgrd = nangrid;
pltgrdt = pltgrd';
ergridt = ergrid';

% Get lat and lon grids
long = linspace(alon1s,alon2s,nxs);
latg = linspace(alat1s,alat2s,nys);

% Mask continents
contmask = NaN(size(pltgrdt));
errmask = NaN(size(ergridt));
lndnan = find(isnan(lndmsk2));
pltgrdt(lndnan) = NaN;
ergridt(lndnan) = NaN;


%% Plot objective map
figure(4)
hold all

pc = pcolor(long,latg,mask(pltgrdt,ergridt,thresh));
set(pc, 'edgecolor','none'); 
plot(coastlon,coastlat,'k','LineWidth',2.5);
plot(bclon,bclat,'k','LineWidth',2.5);
plot(bcslon,bcslat,'k','LineWidth',2.5);
plot(sonlon,sonlat,'k','LineWidth',2.5);
% %% CASE region
xlim([232.0313-360,243.9688-360]);
ylim([30.0438,37.9813]);

LY1 = get(gca,'YLim');
LX1 = get(gca,'XLim');
set(gca,'TickDir','Out')
xticks([-126,-124,-122,-120,-118]);
yticks([31,33,35,37]);
xticklabels(['234';'236';'238';'240';'242']);

caxis([0,2.2]); 
CL = get(gca,'CLim');
cbh = colorbar();
set(cbh,'YTick',linspace(CL(1),CL(2),3)); % Colorbar ticks
set(gca,'dataaspectratio',[1 cos(pi/180*(alat1s+alat2s)/2) 1]);



%% ======= Step 6: Delineate 50% and 80% thresholds ========
pltgrdmsk = mask(pltgrdt,ergridt,thresh); % Mask all non-population values using threshold from obj. maps
zmax = max(pltgrdmsk(:));
grdres = 0.1;

% Threshold: >80%
zpct1 = 80;
zthrs1 = zpct1*0.01*zmax; % Cut-off threshold 
pid1 = pltgrdmsk(find(pltgrdmsk > zthrs1));
[cid1,rid1] = find(pltgrdmsk > zthrs1);
latsmin1 = latg(cid1)-(grdres/2);
latsmax1 = latg(cid1)+(grdres/2);
lonsmin1 = long(rid1)-(grdres/2)+360;
lonsmax1 = long(rid1)+(grdres/2)+360;
thrsh1bnd = boundary(lonsmin1',latsmin1'); % Boundary around all grid cells

% Threshold: >50-80%
zpct2 = 50;
zthrs2 = zpct2*0.01*zmax; % Cut-off threshold 
pid2 = find(pltgrdmsk > zthrs2);
[cid2,rid2] = find(pltgrdmsk > zthrs2);
latsmin2 = latg(cid2)-(grdres/2);
latsmax2 = latg(cid2)+(grdres/2);
lonsmin2 = long(rid2)-(grdres/2)+360;
lonsmax2 = long(rid2)+(grdres/2)+360;
thrsh2bnd = boundary(lonsmin2',latsmin2'); % Boundary around all grid cells
% Gridpoints in Threshold #2 but NOT in Threshold #1
shpunqidx = find(~ismember([lonsmin2',latsmin2'],[lonsmin1',latsmin1'],'rows'));



%% ======= Plot: Grid Cells and Boundary =======
figure(10)

hold on
scatter(lonsmin1,latsmin1,'r');
plot(lonsmin1(thrsh1bnd),latsmin1(thrsh1bnd),'k');
scatter(lonsmin2(shpunqidx),latsmin2(shpunqidx),'b');
plot(lonsmin2(thrsh2bnd),latsmin2(thrsh2bnd),'k');
plot(coastlonall+360,coastlatall,'k','LineWidth',2);
hold off

xlim([232.0313,243.9688]);
ylim([30.0438,37.9813]);

LY1 = get(gca,'YLim');
LX1 = get(gca,'XLim');
set(gca,'TickDir','Out')
xticks([-126+360,-124+360,-122+360,-120+360,-118+360]);
yticks([31,33,35,37]);
xticklabels(['234';'236';'238';'240';'242']);
set(gca,'dataaspectratio',[1 cos(pi/180*(27+40)/2) 1]);

