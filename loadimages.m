% LOADIMAGES  Build an image structure array from image file paths.
%
%   images = loadimages(imglist)
%   images = loadimages(imglist, cam)
%
% Constructs an image structure array and associated camera objects from
% metadata extracted from the listed image file(s). The non-empty fields of the
% provided camera structure (cam) are applied as defaults and scaled as needed
% based on the image size.
%
% Image fields:
% path     - Path to image
% date     - Serial date of image capture
% shutter  - Shutter speed (s)
% aperture - Aperture
% iso      - ISO speed rating
% ev       - Exposure value (https://en.wikipedia.org/wiki/Exposure_value)
%
% GPS fields (if available)
% gps.date       - Serial date of GPS info
% gps.lnglat     - GPS position [longitude latitude]
% gps.altitude   - GPS altitude
% gps.satellites - Number of GPS satellites used
%
% NOTE: Pertinent fields missing in metadata extracted by imfinfo.m:
% camera serial#; lens make, model and serial#; GPS subsec time
%
% Inputs:
%   imglist - List of image paths (string or cell array)
%   cam     - Camera object
%
% Outputs:
%   images  - Image structure array
%
% See also: imfinfo, camera.resize

function images = loadimages(imglist, cam)

  % Check input arguments
  if nargin < 1
    error('Specify list of image paths.');
  end
  if nargin > 1
    if class(cam) ~= 'camera'
      error('Second input must be a camera object.');
    end
  end

  % Convert character string to cell array
  if ischar(imglist)
    imglist = {imglist};
  end

  % Load metadata
  images = struct('path', imglist);
  N = length(imglist);
  for i = 1:N

    % Retrieve metadata
    metadata = imfinfo(imglist{i});

    % Path
    % images(i).fullpath = fullfile(pwd, imglist{i});

    % Date & Time (serial date)
    if isfield(metadata.DigitalCamera, 'DateTimeOriginal')
      DateTimeOriginal = deblank(metadata.DigitalCamera.DateTimeOriginal);
      if isfield(metadata.DigitalCamera, 'SubsecTimeOriginal')
        SubsecTime = deblank(strtrim(metadata.DigitalCamera.SubsecTimeOriginal));
      elseif isfield(metadata.DigitalCamera, 'SubsecTime')
        SubsecTime = deblank(strtrim(metadata.DigitalCamera.SubsecTime));
      else
        SubsecTime = '00';
      end
      % images(i).date = [DateTimeOriginal '.' SubsecTime];
      images(i).date = datenum([DateTimeOriginal '.' SubsecTime], 'yyyy:mm:dd HH:MM:SS.FFF');
    end

    % Exposure
    % Shutter speed (exposure time in seconds)
    if isfield(metadata.DigitalCamera, 'ExposureTime')
      images(i).shutter = metadata.DigitalCamera.ExposureTime; % (in seconds)
    end
    % Aperture (FNumber)
    if isfield(metadata.DigitalCamera, 'FNumber')
      images(i).aperture = metadata.DigitalCamera.FNumber;
    end
    % ISO
    if isfield(metadata.DigitalCamera, 'ISOSpeedRatings')
      images(i).iso = metadata.DigitalCamera.ISOSpeedRatings;
    end
    % Exposure value
    % Logarithmic base-2 scale. EV 0 corresponds to an exposure time of 1
    % second with an aperture of 1.0 at ISO 100.
    if isfield(metadata.DigitalCamera, {'ExposureTime', 'FNumber', 'ISOSpeedRatings'})
      images(i).ev = log2(100 * images(i).aperture^2 / (images(i).iso * images(i).shutter));
    end

    % GPS
    if isfield(metadata, 'GPSInfo')

      % GPS Date & Time (serial date)
      if isfield(metadata.GPSInfo, 'GPSDateStamp')
        GPSDateStamp = deblank(strtrim(metadata.GPSInfo.GPSDateStamp));
        HMS = metadata.GPSInfo.GPSTimeStamp;
        gpsdatestr = [GPSDateStamp ' ' sprintf('%02d', HMS(1)) ':' sprintf('%02d', HMS(2)) ':' sprintf('%02d', HMS(3))];
        % images(i).gpsdate = gpsdatestr;
        images(i).gps.date = datenum(gpsdatestr, 'yyyy:mm:dd HH:MM:SS');
        % NOTE: Missing GPS SubsecTime field
      end

      % GPS Position (lng, lat, altitude)
      if isfield(metadata.GPSInfo, {'GPSLatitude','GPSLongitude'})
        % Latitude (decimal degrees)
        lat = metadata.GPSInfo.GPSLatitude;
        lat = lat(1) + lat(2) / 60 + lat(3) / 3600;
        if deblank(strtrim(metadata.GPSInfo.GPSLatitudeRef)) == 'S'
          lat = -lat;
        end
        % Longitude (decimal degrees)
        lng = metadata.GPSInfo.GPSLongitude;
        lng = lng(1) + lng(2) / 60 + lng(3) / 3600;
        if deblank(strtrim(metadata.GPSInfo.GPSLongitudeRef)) == 'W'
          lng = -lng;
        end
        images(i).gps.lnglat = [lng lat];
        % if exist('wgs2utm',' file')
        %   [xutm, yutm, utmzone] = wgs2utm(latitude,longitude);
        %   images(i).gps_utm = [xutm yutm utmzone];
        % end
      end

      % GPS Altitude (FIXME: meters? relative to what surface?)
      if isfield(metadata.GPSInfo, 'GPSAltitude')
        images(i).gps.altitude = metadata.GPSInfo.GPSAltitude;
      end

      % GPS Satellites (number of satellites)
      if isfield(metadata.GPSInfo, 'GPSSatellites')
        images(i).gps.satellites = str2double(deblank(strtrim(metadata.GPSInfo.GPSSatellites)));
      end
    end

    % Camera parameters
    % Image size in pixels [ny, nx]
    if nargin > 1
      if isempty(cam.imgsz)
        imgcam = cam;
        imgcam.imgsz = [metadata.Width metadata.Height];
      else
        imgcam = cam.resize([metadata.Width metadata.Height]);
      end
    else
      imgcam = camera('imgsz', [metadata.Width metadata.Height]);
    end
    % Sensor size [mm width, mm height]
    if nargin > 1 && isempty(imgcam.sensorsz)
      % Look up using camera make and model
      if isfield(metadata, {'Make', 'Model'})
        imgcam.sensorsz = sensorSize(metadata.Make, metadata.Model);
        % images(i).makemodel = [deblank(strtrim(metadata.Make)) ' ' deblank(strtrim(metadata.Model))]
      end
    end
    if ~isempty(imgcam.sensorsz) && isempty(imgcam.f)
      % Focal length in mm
      if isfield(metadata.DigitalCamera, 'FocalLength')
        imgcam.fmm = metadata.DigitalCamera.FocalLength;
      elseif isfield(metadata.DigitalCamera, 'FocalLengthIn35mmFilm')
        imgcam.fmm = metadata.DigitalCamera.FocalLengthIn35mmFilm;
        % FIXME: Convert to true focal length using sensor size?
        warning('True focal length not found, using 35mm equivalent.');
      end
    end
    images(i).cam = imgcam;

  end

  %% BACKUP

  % % Sort by capture time
  % [~, order] = sort([images(:).datenum]);
  % images = images(order);
  %
  % % Check existence of non-default metadata fields
  % fields = {'Make', 'Model'};
  % missing = find(~isfield(metadata, fields));
  % for m = 1:length(missing)
  %   metadata.(fields{m}) = '';
  % end
