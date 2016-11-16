function horizon = sunMask( location, dem, varargin )
%Retrieve local horizon obstruction from Digital Elevation Model (DEM)
%
%   Summary:
%       Compute local horizon due to terrain surroundings and topography
%       using a Digital Elevation Model tile. By default, this algorithm
%       uses the SRTM-based DEM definition (WGS84 datum), whose tiles (with
%       a 3 arc-second resolution) can be downloaded at :
%       http://srtm.csi.cgiar.org/SELECTION/inputCoord.asp. However, users
%       can apply this algorithm to any DEM whose reference ellipsoid is
%       already implemented (see 'getEllipsoidProperties' function) or can
%       define a system of their own (basically by defining semi major axis
%       and flattening of the given ellipsoid). The only requirement is
%       that the DEM elevation data must be defined in the geographic
%       coordinate system.
%       
%       Horizon is defined as a vector of elevation (or zenith) points (in
%       °) with their corresponding azimuth (in °). Typically, it can then
%       be used to assess topography effects into solar radiation
%       computation ("sun mask"). The present algorithm has been depicted
%       in the following article: "Pillot, B., Muselli, M., Poggi, P.,
%       Haurant, P., Dias, J.B., Development and validation of a new
%       efficient SRTM DEM-based horizon model combined with optimization
%       and error prediction methods, Solar Energy, 2016, vol. 129, pp.
%       101-115." It is both faster and more accurate than other typical
%       DEM-derived horizon models, such as Carnaval or r.horizon (plugin
%       of the GRASS GIS software). If you intend to publish material using
%       this algorithm, we kindly ask you to cite this work.
%
%   Syntax:
%       * [elevation, azimuth] = SUNMASK(location, dem)
%       * [elevation, azimuth] = SUNMASK(location, dem, distance)
%       * [elevation, azimuth] = SUNMASK(location, dem, distance, precision)
%       * [elevation, azimuth] = SUNMASK(location, dem, distance, precision, ellipsoidModel)
%
%   Description:
%       * SUNMASK(location, dem) - compute horizon, i.e. angular elevation,
%       zenith and azimuth of surrounding topography points, according to
%       the given location and geo-referenced Digital Elevation Model (DEM)
%       * SUNMASK(location, dem, distance) - compute horizon with respect
%       to viewing distance defined by user (in km)
%       * SUNMASK(location, dem, distance, precision) - compute horizon
%       with specific precision, i.e. azimuth step between points composing
%       it
%       * SUNMASK(location, dem, distance, precision, ellipsoidModel) -
%       compute horizon with respect to the given ellipsoid model used to
%       represent DEM elevation data
%    
%   Inputs:
%       * LOCATION (required) - structure with both latitude and longitude
%       fields corresponding to latitude and longitude of the given
%       location:
%           * LOCATION.LATITUDE - (scalar) geographic latitude of the given
%           location (in °)
%           * LOCATION.LONGITUDE - (scalar) geographic longitude of the
%           given location (in °)
%       * DEM (required) - structure with 3 fields corresponding to
%       elevation, latitude and longitude of the DEM points:
%           * DEM.ELEVATION - (matrix) SRTM-based Digital Elevation Model¹
%           tile defined as a matrix array of elevation points (in m) (can
%           be downloaded as a GeoTiff image at
%           http://srtm.csi.cgiar.org/SELECTION/inputCoord.asp).
%           * DEM.LATITUDE - (vector) geographic latitude of the DEM points
%           (in °), with number of elements equal to size(DEM.ELEVATION, 1)
%           * DEM.LONGITUDE - (vector) geographic longitude of the DEM
%           points (in °), with number of elements eqal to
%           size(DEM.ELEVATION, 2)
%       * DISTANCE (optional) - (positive scalar) horizon viewing distance
%       (in km); default = 50 km
%       * PRECISION (optional) - (positive scalar) azimuth precision of the
%       computing horizon (in °); default = 1°
%       * ELLIPSOIDMODEL (optional) - corresponding ellipsoid model used to
%       represent the DEM elevation data; default = 'WGS84' (SRTM DEM
%       definition)
%
%       ¹Actually, in the case of the SRTM-based DEM, the altitude baseline
%       is the EGM96 geoid. Though undulations of the EGM96 with respect to
%       the WGS84 ellipsoid remain low (between -107m and 85m), the user
%       can take them into consideration by downloading the NGA model at:
%       http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/egm96.html
%       In the case of another reference ellipsoid, be sure elevation data
%       are given as ellipsoidal heights, or use the corresponding gravity
%       model for taking into account the geoid undulation.
%
%   Outputs:
%       * HORIZON - structure with 3 fields corresponding to elevation,
%       zenith and azimuth of the horizon points:
%           * HORIZON.ELEVATION - (vector) elevation angle (in °) of the
%           points composing the surrounding horizon
%           * HORIZON.ZENITH - (vector) zenith angle (in °) of the points
%           composing the surrounding horizon
%           * HORIZON.AZIMUTH - (vector) azimuth (in °) of the points
%           composing the horizon
%
%   Examples:
%       * See the script 'Example/horizonExample.m' for how to use this
%       function
%
%   Other m-files required:
%       * GETELLIPSOIDPROPERTIES.m
%
%   Subfunctions:
%       * GETELEVATIONANGLE
%       * GETAZIMUTH
%       * GETISOMETRICLATITUDE
%       * GEOGRAPHIC2CARTESIAN
%       * DIST2DEG
%       * GEODESIC
%       * LOCALSPHERERADIUS
%
%   See also:
%       * GETELLIPSOIDPROPERTIES

%   Author: Dr. Benjamin Pillot
%   Address: Universidade do Vale do Rio dos Sinos (Unisinos), São 
%   Leopoldo, RS, Brazil
%   email: benjaminfp@unisinos.br
%   Website: http://www.
%   Date: 02-Oct-2016; Last revision: 02-Oct-2016
%
%   Copyright (c) 2016, Benjamin Pillot
%   All rights reserved.


%% Check input arguments
narginchk(2, 5);

% Default horizon viewing distance (50 km), precision (1°) and ellipsoid
% model (WGS84)
distance = 50;
precision = 1;
ellipsoidModel = 'WGS84';

if nargin == 3
    distance = varargin{1};
end
if nargin == 4
    precision = varargin{2};
end
if nargin == 5
    ellipsoidModel = varargin{3};
end

assert(isnumeric(distance) && distance > 0 && distance <= 500 && isscalar(distance), ...
    'distance must be a positive scalar <= 500 km');
assert(isnumeric(precision) && precision > 0 && isscalar(precision), ...
    'precision must be a positive scalar');
assert(ischar(ellipsoidModel) && isrow(ellipsoidModel), 'ellipsoid model must be a string');


%% Initialization
% DEM resolution
demResolution = abs(dem.latitude(1) - dem.latitude(2));

% Retrieve the DEM pixel corresponding to location
[~, yPosition] = min(abs(dem.latitude - location.latitude));
[~, xPosition] = min(abs(dem.longitude - location.longitude));

% Maximum computation distance
[latitudeDistance, longitudeDistance] = dist2deg(distance, location, ellipsoidModel);
topoLatDistance = round(latitudeDistance/demResolution);      % Distance between observer and the further topography point located on a latitude isoline
topoLonDistance = round(longitudeDistance/demResolution);     % Distance between observer and the further topography point located on a longitude isoline


% Take into consideration possible out mapping : if maximum distance is out
% of the dem area, limit to the dem area
northAreaDistance = topoLatDistance - abs(min(0, yPosition - topoLatDistance - 1));
southAreaDistance = topoLatDistance - max(0, yPosition + topoLatDistance - size(dem.elevation, 1));
westAreaDistance = topoLonDistance - abs(min(0, xPosition - topoLonDistance - 1));
eastAreaDistance = topoLonDistance - max(0, xPosition + topoLonDistance - size(dem.elevation, 2));    

% Retrieve specific dem area for horizon computation
studyArea = dem.elevation(yPosition - northAreaDistance : yPosition + southAreaDistance, ...
    xPosition - westAreaDistance : xPosition + eastAreaDistance);
[ySize, xSize] = size(studyArea);

% Geo-referencing meshgrid
[longitude, latitude] = meshgrid(dem.longitude(xPosition - (westAreaDistance + 1) : xPosition + eastAreaDistance), ...
    dem.latitude(yPosition - (northAreaDistance + 1) : yPosition + southAreaDistance));


% Azimuth of all pixel corners seen from observer point
xObs = westAreaDistance + 1;
yObs = northAreaDistance + 1;

azimuth = (180/pi) * getAzimuth((pi/180)*location.latitude, (pi/180)*location.longitude, ...
    (pi/180)*(latitude - demResolution/2), (pi/180)*(longitude + demResolution/2), ellipsoidModel);

% Specific azimuth values for NE (-180 to 90) and NW (90 to 180) areas
azimuthNE = azimuth;
azimuthNE(1 : yObs, xObs) = azimuthNE(1 : yObs, xObs) - 360;

azimuthNW = azimuth;
azimuthNW(1 : yObs, xObs + 1) = azimuthNW(1 : yObs, xObs + 1) + 360;

% Corresponding elevation angle
zObs = studyArea(yObs, xObs);
elevation = getElevationAngle(zObs, studyArea, (pi/180)*location.latitude, ...
    (pi/180)*latitude(2 : end, 2 : end), (pi/180)*location.longitude, ...
    (pi/180) * longitude(2 : end, 2 : end), ellipsoidModel);

% Elevation vector length
lengthElevation = (90 + precision)/precision;

% Initialization
elevationNE = zeros(northAreaDistance, lengthElevation);
elevationE = zeros(eastAreaDistance, 2*lengthElevation - 1);
elevationS = zeros(southAreaDistance, 2*lengthElevation - 1);
elevationW = zeros(westAreaDistance, 2*lengthElevation - 1);
elevationNW = zeros(northAreaDistance, lengthElevation);

azNE = -180 : precision : -90;
azE = -180 : precision : 0;
azS = -90 : precision : 90;
azW = 0 : precision : 180;
azNW = 90 : precision : 180;


%% Main computation
% Retrieve all elevation angles for iso-azimuth lines (loxodromes)
for isoline = 1 : northAreaDistance
        
    [~, k] = histc(azNE, azimuthNE(isoline + 1, xObs : xSize + 1));
    elevationNE(isoline, k ~= 0) = elevation(isoline, xObs - 1 + k(k ~= 0));
    
    
    [~, k2] = histc(azNW, azimuthNW(isoline + 1, 1 : xObs + 1));
    elevationNW(isoline, k2 ~= 0) = elevation(isoline, k2(k2 ~= 0));
    
end

for isoline = 1 : southAreaDistance
    
    [~, k] = histc(azS, azimuth(yObs + isoline, xSize + 1 : -1 : 1));
    elevationS(isoline, k ~= 0) = elevation(yObs + isoline, end + 1 - k(k ~= 0));
    
end

for isoline = 1 : eastAreaDistance
    
    [~, k] = histc(azE, azimuth(1 : ySize + 1, xObs + isoline));
    elevationE(isoline, k ~= 0) = elevation(k(k ~= 0), xObs + isoline);
    
end

for isoline = 1 : westAreaDistance
    
    [~, k] = histc(azW, azimuth(ySize + 1 : -1 : 1, isoline + 1));
    elevationW(isoline, k ~= 0) = elevation(end + 1 - k(k ~= 0), isoline);
    
end

sunMaskNE = max(elevationNE);
sunMaskE = max(elevationE);
sunMaskS = max(elevationS);
sunMaskW = max(elevationW);
sunMaskNW = max(elevationNW);

% Global azimuth (North to North) and sun mask (elevation angle)
azNtoN = [azNE azE azS azW azNW];
sunMask = [sunMaskNE sunMaskE sunMaskS sunMaskW sunMaskNW];

horizon.elevation = zeros(1, (360 + precision)/precision);

% Corresponding azimuth (from -180° to 180°)
az = -180 : precision : 180;

for r = 1 : length(az)
                
    horizon.elevation(r) = max(sunMask(azNtoN == az(r)));
                            
end


%% Finalization
% Set horizon
horizon.elevation(horizon.elevation < 0) = 0;
horizon.zenith = 90 - horizon.elevation;
horizon.azimuth = 0 : precision : 360;  % Store azimuth from 0 to 360 (for using with solar azimuth functions)

end




%% Sub functions
%   * GETELEVATIONANGLE - compute angular elevation between 2 points of any
%   ellipsoid using ellipsoidal height and geographic coordinates (latitude
%   and longitude)
%   * GETAZIMUTH - compute azimuth between 2 points of the given ellipsoid
%   * GETISOMETRICLATITUDE - compute isometric latitude from geographic
%   latitude
%   * GEOGRAPHIC2CARTESIAN - convert geographic coordinates to cartesian
%   coordinates, with respect to the ellipsoid model
%   * DIST2DEG - convert metric distance into angular distance using local
%   sphere approximation
%       * GEODESIC - compute geodesic on ellipsoid using local sphere
%       approximation
%       * LOCALSPHERERADIUS - retrieve radius of the local sphere
%       approximation, tangent to the ellipsoid


function alpha = getElevationAngle(h_A, h_B, latitude_A, latitude_B, longitude_A, longitude_B, ellipsoidModel)

% Compute cartesian coordinates of point A and B located at altitude
% h_A and h_B from the ellipsoid surface (ellipsoidal heights)
[x_A, y_A, z_A] = geographic2cartesian(latitude_A, longitude_A, h_A, ellipsoidModel);
[x_B, y_B, z_B] = geographic2cartesian(latitude_B, longitude_B, h_B, ellipsoidModel);

% Scalar product between AB and normal to the point A
innerProduct = (x_B - x_A) .* cos(longitude_A) .* cos(latitude_A) + (y_B - y_A) .* ...
    sin(longitude_A) .* cos(latitude_A) + (z_B - z_A) .* sin(latitude_A);

% Angular elevation computation
norm = sqrt((x_B - x_A).^2 + (y_B - y_A).^2 + (z_B - z_A).^2);
alpha = asin(innerProduct./norm) * 180/pi;

end


function az = getAzimuth(lat1, lon1, lat2, lon2, ellipsoidModel)

% Retrieve isometric latitudes
L1 = getIsometricLatitude(lat1, ellipsoidModel);
L2 = getIsometricLatitude(lat2, ellipsoidModel);

% Compute azimuth
az = atan2((lon1 - lon2), (L1 - L2));        
        
end


function isometricLatitude = getIsometricLatitude(latitude, ellipsoidModel)

[~, ~, ~, e] = getEllipsoidProperties(ellipsoidModel);
term1 = tan((pi/4) + (latitude/2));
num = 1 - e * sin(latitude);
denom = 1 + e * sin(latitude);
term2 = (num./denom).^(e/2);

% Result
isometricLatitude = log(term1 .* term2);
        
end


function [x, y, z] = geographic2cartesian(latitude, longitude, altitude, ellipsoidModel)

% Compute semi-major axis and eccentricity of the specified ellipsoid
[a, ~, ~, e] = getEllipsoidProperties(ellipsoidModel);

% Compute ellipsoid normal
N = a./sqrt(1 - (e^2)*(sin(latitude)).^2);

% Compute cartesian coordinates from geographic coordinates
x = (N + altitude) .* cos(longitude) .* cos(latitude);
y = (N + altitude) .* sin(longitude) .* cos(latitude);
z = (N * (1 - e^2) + altitude) .* sin(latitude);        
        
end


function [deltaLat, deltaLon] = dist2deg(distance, location, ellipsoidModel)

% Initialization
distanceEps = 1e-2;
latMin = 0;
latMax = 10;
lonMin = 0;
lonMax = 10;

% Computation (dichotomy)
while true
        
    deltaLat = (latMin + latMax)/2;
    deltaLon = (lonMin + lonMax)/2;
    
    distVarLat = geodesic(pi/180*location.latitude, pi/180*(location.latitude + deltaLat), 0, ellipsoidModel);
    distVarLon = geodesic(pi/180*location.latitude, pi/180*location.latitude, pi/180*deltaLon, ellipsoidModel);
    
    if abs(distVarLat - distance) < distanceEps && abs(distVarLon - distance) < distanceEps
        break       
    end
    
    if distVarLat < distance        
        latMin = deltaLat;        
    else        
        latMax = deltaLat;        
    end
    
    if distVarLon < distance        
        lonMin = deltaLon;        
    else        
        lonMax = deltaLon;        
    end    
    
end

end


function arc_AB = geodesic(lat1, lat2, deltaLon, ellipsoidModel)

% Compute local sphere radius
R = localSphereRadius(lat1, ellipsoidModel);

% Compute angle
gamma = acos(cos(lat1).*cos(lat2).*cos(deltaLon) + sin(lat1).*sin(lat2));

% Compute great circle
arc_AB = R .* gamma;

end


function R = localSphereRadius(latitude, ellipsoidModel)

% Get ellipsoid definition with respect to given model
[a, ~, ~, e] = getEllipsoidProperties(ellipsoidModel);

p = ((a * (1 - e^2))./((1 - (e^2)*(sin(latitude)).^2).^(3/2))) * 10^-3;
N = (a./sqrt(1 - (e^2) * (sin(latitude)).^2)) * 10^-3;

% Compute radius
R = sqrt(N.*p);

end
