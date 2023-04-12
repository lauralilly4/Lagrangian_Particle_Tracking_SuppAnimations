%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Code for Supplemental Animations to:
%%%     Lilly et al. Using a Lagrangian particle tracking model to evaluate 
%%%     impacts of El Niño-related advection on euphausiids in the southern 
%%%     California Current System
%%% Author: Laura E. Lilly
%%% Updated: 27 Jan 2022

% CODING STEP 1: Load and read files: US/Baja California coastline and land mask


%% ======= Load files: coastline and land mask =======
% Coastline shapefiles
ca_coast = shaperead('Western.shp');
bc_coast = shaperead('bc_entidad.shp');
bcs_coast = shaperead('bcs_entidad.shp');
son_coast = shaperead('son_entidad.shp');
% Land mask
[X,Y,Z] = grdread2('ETOPO1_Ice_g_gmt4.grd');
Zd = double(Z); % Convert grid values to class 'double' (from int32)


%% ======= Read in coastline shapefiles =======
% California
coastlat = [];
coastlon = [];

calat = ca_coast.Y;
for c = 1:50000
    clat = ca_coast(c).Y;
    clon = ca_coast(c).X;
    coastlat = [coastlat, clat];
    coastlon = [coastlon, clon];
end

% Baja California
bclat = bc_coast.Y;
bclon = bc_coast.X;

bclat(find(bclat < 28)) = NaN;
bclon(find(bclat < 28)) = NaN;

bclat(find(bclat > 32.5)) = NaN;
bclon(find(bclat > 32.5)) = NaN;

bclat(find(bclat > 31.7 & bclon > -115.2)) = NaN;
bclon(find(bclat > 31.7 & bclon > -115.2)) = NaN;

% Baja California Sur
bcslat = bcs_coast.Y;
bcslon = bcs_coast.X;

bcslat(find(bcslat > 27.99)) = NaN;
bcslon(find(bcslat > 27.99)) = NaN;

% Sonora
sonlat = son_coast.Y;
sonlon = son_coast.X;
sonlat(find(sonlat > 31.8)) = NaN;
sonlon(find(sonlat > 31.8)) = NaN;

sonlat(find(sonlat > 31.5 & sonlon > -113)) = NaN;
sonlon(find(sonlat > 31.5 & sonlon > -113)) = NaN;

% Combine all sections
coastlatall = [coastlat,bclat,bcslat,sonlat];
coastlonall = [coastlon,bclon,bcslon,sonlon];


%% ======= Read in files to create land mask =======
lonmn = -126.717; % 'alon1s' from CalCOFI obj mapping code! Use same value to convert these!!
lonmx = -112.117;
latmn = 24.15; % 'alat1s' - same as above
latmx = 37.967; %'alat2s' 

[xvls,xids] = find(X >= lonmn & X <= lonmx);
[yvls,yids] = find(Y >= latmn & Y <= latmx);

xcut = X(xids);
ycut = Y(yids);
zcut = Zd(yids,xids);

% Convert distances to km
tokm = 111.324; % Lat -> km 
tokmt = tokm*cos(pi/180*(latmn+latmx)/2); % Lon -> km
exs = tokmt*(xcut-lonmn); % Convert 'etopo' values to km
eys = tokm*(ycut-latmn);

xrs = tokmt*(lonmx-lonmn); % Conversion of total x distance to km
yrs = tokm*(latmx-latmn);
nxs = 100;
nys = 101;
xgbins = linspace(0,xrs,nxs); % Vector of gridpoints: min (0), max (diff), n steps between
ygbins = linspace(0,yrs,nys); 

lndmsk = NaN(nys,nxs);

for zx = 1:nxs
    for zy = 1:nys
        dxg = (xgbins(2)-xgbins(1))/2;
        dyg = (ygbins(2)-ygbins(1))/2;
        [~,xi] = find( exs > (xgbins(zx)-dxg) & exs <= (xgbins(zx)+dxg) );
        [~,yi] = find( eys > (ygbins(zy)-dyg) & eys <= (ygbins(zy)+dyg) );

        zxy = zcut(yi,xi);
        zmed = median(zxy(:));
        
        if zmed > 0
            medin = NaN;
        else
            medin = zmed;
        end
        lndmsk(zy,zx) = medin;
    end
end
   

% Cut-off lat & lon for Salton Sea and Gulf of California
gclx = tokmt*(-114-lonmn); % Gulf of California
gcly = tokm*(27.5-latmn);
gcux = tokmt*(-115-lonmn);
gcuy = tokm*(29.5-latmn); % GoCA - upper
ssx = tokmt*(-117-lonmn);
ssy = tokm*(32-latmn);
        
lndmsk2 = NaN(nys,nxs);

for lx = 1:nxs
    for ly = 1:nys
        if xgbins(lx) > gclx && ygbins(ly) > gcly
            lndmsk2(ly,lx) = NaN;
        elseif xgbins(lx) > gcux && ygbins(ly) > gcuy
            lndmsk2(ly,lx) = NaN;
        elseif xgbins(lx) > ssx && ygbins(ly) > ssy
            lndmsk2(ly,lx) = NaN;
        else
            lndmsk2(ly,lx) = lndmsk(ly,lx);
        end
    end
end


