function [semiMajorAxis, flattening, semiMinorAxis, eccentricity] = getEllipsoidProperties( model )
%Retrieve ellipsoid properties from model
%
%   Summary:
%       Retrieve ellipsoid parameters from model. At present time, only
%       WGS84, ETRS89 and IERS are defined, but user can easily add other
%       ellipsoid definitions (semi major axis and flattening)
%
%   Syntax:
%       * [semiMajorAxis, flattening, semiMinorAxis, eccentricity] = ELLIPSOIDMODEL(model)
%
%   Description:
%       * ELLIPSOIDMODEL(model) - return specific ellipsoid parameters
%    
%   Inputs:
%       * MODEL - string corresponding to the model¹ name ('WGS84', etc.)
%
%       ¹Currently implemented models:
%       - WGS84, ETRS89, IERS
%
%   Outputs:
%       * SEMIMAJORAXIS - semi major axis of the given ellipsoid
%       * FLATTENING - flattening of the given ellipsoid
%       * SEMIMINORAXIS - sem minor axis of the given ellipsoid
%       * ECCENTRICITY - eccentricity of the given ellipsoid
%
%   Examples:
%       * getEllipsoidProperties('WGS84') retrieves all fundamental
%       parameters (semi-major and semi-minos axes, flattening and
%       eccentricity) of the corresponding reference ellipsoid World
%       Geodetic System 1984 (NGA, 2000. Department of Defense World
%       Geodetic System 1984 – its definition and relationships with local
%       geodetic systems. Tech. Rep., National Geospatial-Intelligence
%       Agency, Springfield, USA.)
%
%   See also:
%       * SUNMASK

%   Author: Dr. Benjamin Pillot
%   Address: Universidade do Vale do Rio dos Sinos (Unisinos), São 
%   Leopoldo, RS, Brazil
%   email: benjaminfp@unisinos.br
%   Website: http://www.
%   Date: 03-Oct-2016; Last revision: 03-Oct-2016
%
%   Copyright (c) 2016, Benjamin Pillot
%   All rights reserved.

% Check input arguments
narginchk(0, 1);

if nargin == 0    
    model = 'WGS84';
end

model = validatestring(model, {'WGS84', 'ETRS89', 'IERS'});
% model = validatestring(model, {'WGS84', 'ETRS89', 'IERS', 'your_model'});

switch model
    
    case 'WGS84'    % World Geodetic System 1984
        
        flattening = 1/298.257223563;
        semiMajorAxis = 6378137;
        
    case 'ETRS89'   % European Terrestrial Reference System 1989
        
        flattening = 1/298.257222101;
        semiMajorAxis = 6378137;
        
    case 'IERS'     % International Earth Rotation (2003)
        
        flattening = 1/298.25642;
        semiMajorAxis = 6378136.6;
        
    %case 'your_model'
    %   Add here parameters of your ellipsoid model (flattening and semi
    %   major axis)
        
end

% Computation of semi-minor axis and eccentricity
semiMinorAxis = semiMajorAxis - semiMajorAxis * flattening;
eccentricity = sqrt(1 - (semiMinorAxis^2/semiMajorAxis^2));


end
