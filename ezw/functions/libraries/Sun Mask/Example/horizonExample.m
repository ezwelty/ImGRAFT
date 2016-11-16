%HORIZONEXAMPLE - Example of horizon computation using a 3'' SRTM DEM
%
%   Summary:
%   	Compute horizon using 3'' SRTM DEM (45°N-40°N; 5°E-10°E) downloaded
%   	at http://srtm.csi.cgiar.org/SELECTION/inputCoord.asp
%       Return horizon computed for 3 different conditions (default, 100 km
%       viewing distance, 0.5° precision) as well as the corresponding
%       computation times.
%
%   Description:
%       * HORIZONEXAMPLE - return 3 horizon structures for location 42°N,
%       9°E and their corresponding computation time
%
%   Outputs:
%       * HORIZON, COMPUTATIONTIME
%       * HORIZON_DISTANCE_100, COMPUTATIONTIME_DISTANCE_100
%       * HORIZON_PRECISION_0_5, COMPUTATIONTIME_PRECISION_0_5
%
%   MAT-files required:
%       * INPUTS.MAT
%
%   See also:
%	* SUNMASK, GETELLIPSOIDPROPERTIES

%   Author: Dr. Benjamin Pillot
%   Address: Universidade do Vale do Rio dos Sinos (Unisinos), São 
%   Leopoldo, RS, Brazil
%   email: benjaminfp@unisinos.br
%   Website: http://www.
%   Date: 03-Oct-2016; Last revision: 03-Oct-2016
%
%   Copyright (c) 2016, Benjamin Pillot

% Load DEM and location
load inputs.mat

% 1. Default function call
time1 = tic;
horizon = sunMask(location, dem);
computationTime = toc(time1);

% 2. Call function with different viewing distance (100km)
time2 = tic;
horizon_distance_100 = sunMask(location, dem, 100);
computationTime_distance_100 = toc(time2);

% 3. Call function with different precision (0.5°)
time3 = tic;
horizon_precision_0_5 = sunMask(location, dem, 50, 0.5);
computationTime_precision_0_5 = toc(time3);

clear time1 time2 time3