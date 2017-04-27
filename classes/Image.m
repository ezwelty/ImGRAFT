classdef Image < handle
% IMAGE Photographic image data structure.
%
% Image Properties:
% cam      - Camera object
%
% Image Properties (read-only):
% info     - File information from imfinfo
% file     - Path to image file
% date_str - Capture date and time as a string ('yyyy-mm-dd HH:MM:SS.FFF')
% date_num - Cature date and time as a serial date number
% shutter  - Shutter speed in seconds
% aperture - Lens aperture
% iso      - Film speed
% ev       - Exposure value
% gps      - GPS metadata
% size     - Size of original image [nx|ncols|width, ny|nrows|height]
%
% Image Properties (dependent):
% scale    - Scaling between original and camera image size
%
% Image Methods:
% Image - Construct a new Image object
%
% Image Methods (cached):
% read  - Read image data from file
% clear - Clear object cache
%
% See also imfinfo, datestr, datenum

  properties
    cam
  end

  properties (SetAccess = private)
    info
  end

  properties (Dependent, SetAccess = private)
    file
    date_str
    date_num
    shutter
    aperture
    iso
    ev
    gps
    size
  end
    
  properties (Dependent)
    scale
  end
  
  properties (SetAccess = private, Hidden = true)
    cached_I
  end
  
  methods

    % Image creation
    
    function images = Image(files, cam)
    % IMAGE  Construct a new Image object.
    %
    %   img = Image(files, cam = Camera())
    %
    % Image size, sensor size, and focal length are loaded from the file
    % unless overloaded by cam.
    %
    % Inputs:
    %   file - Path to image file
    %   cam  - Camera object

      % Check inputs
      if nargin < 1
        return
      end
      if isempty(files)
        temp = Image();
        images = temp(false);
        return
      end
      if ~isa(files, 'cell')
        files = {files};
      end
      if nargin < 2
        cam = Camera();
      end
      if nargin > 1 && ~isa(cam, 'Camera')
        error('Not an object of class Camera.');
      end
      % Expand paths
      files = cellfun(@expand_path, files, 'UniformOutput', false);
      files = [files{:}];
      % Preallocate array
      images(length(files)) = Image();
      for i = 1:length(files)
        images(i) = foreach(images(i), files{i}, cam);
      end
      % For each file
      function img = foreach(img, file, cam)
        % Metadata
        img.info = imfinfo(file);
        % Camera
        img.cam = cam;
        % Image size in pixels [nx, ny]
        if isempty(img.cam.imgsz)
          img.cam.imgsz = [img.info.Width, img.info.Height];
        else
          % Check that target size is compatible with image size
          Camera.getScaleFromSize([img.info.Width, img.info.Height], img.cam.imgsz);
        end
        % Sensor size [mm width, mm height]
        if isempty(img.cam.sensorsz) && all(isfield(img.info, {'Make', 'Model'}))
          img.cam.sensorsz = Camera.sensorSize(img.info.Make, img.info.Model);
        end
        % Focal length in mm (if not already set in pixels)
        if ~isempty(img.cam.sensorsz) && isempty(img.cam.f)
          if isfield(img.info.DigitalCamera, 'FocalLength')
            img.cam.fmm = img.info.DigitalCamera.FocalLength;
          elseif isfield(img.info.DigitalCamera, 'FocalLengthIn35mmFilm')
            % FIXME: Convert to true focal length using sensor size?
            img.cam.fmm = img.info.DigitalCamera.FocalLengthIn35mmFilm;
            warning('True focal length not found, using 35mm equivalent.');
          end
        end
      end
    end

    function value = get.file(img)
      value = img.info.Filename;
    end

    function value = get.date_str(img)
      if isfield(img.info.DigitalCamera, 'DateTimeOriginal')
        date_time = strsplit(deblank(img.info.DigitalCamera.DateTimeOriginal), ' ');
        date_time{1} = strrep(date_time{1}, ':', '-');
        value = [date_time{1} ' ' date_time{2}];
        if isfield(img.info.DigitalCamera, 'SubsecTimeOriginal')
          value = [value '.' deblank(strtrim(img.info.DigitalCamera.SubsecTimeOriginal))];
        elseif isfield(img.info.DigitalCamera, 'SubsecTime')
          value = [value '.' deblank(strtrim(img.info.DigitalCamera.SubsecTime))];
        end
      else
        value = [];
      end
    end

    function value = get.date_num(img)
      format = 'yyyy-mm-dd HH:MM:SS';
      if ~isempty(strfind(img.date_str, '.'))
        format = [format '.FFF'];
      end
      value = datenum(img.date_str, format);
    end

    function value = get.shutter(img)
      if isfield(img.info.DigitalCamera, 'ExposureTime')
        value = img.info.DigitalCamera.ExposureTime; % (in seconds)
      else
        value = [];
      end
    end

    function value = get.aperture(img)
      if isfield(img.info.DigitalCamera, 'FNumber')
        value = img.info.DigitalCamera.FNumber;
      else
        value = [];
      end
    end

    function value = get.iso(img)
      if isfield(img.info.DigitalCamera, 'ISOSpeedRatings')
        value = img.info.DigitalCamera.ISOSpeedRatings;
      else
        value = [];
      end
    end

    function value = get.ev(img)
      % https://en.wikipedia.org/wiki/Exposure_value
      value = log2(1000 * img.aperture^2 / (img.iso * img.shutter));
    end

    function value = get.gps(img)
      value = [];
      if isfield(img.info, 'GPSInfo')
        % GPS Date & Time (serial date)
        % NOTE: Missing GPS SubsecTime field.
        if isfield(img.info.GPSInfo, 'GPSDateStamp')
          GPSDateStamp = deblank(strtrim(img.info.GPSInfo.GPSDateStamp));
          HMS = img.info.GPSInfo.GPSTimeStamp;
          value.datestr = [GPSDateStamp ' ' sprintf('%02d', HMS(1)) ':' sprintf('%02d', HMS(2)) ':' sprintf('%02d', HMS(3))];
          value.datenum = datenum(value.datestr, 'yyyy:mm:dd HH:MM:SS');
        end
        % GPS Position (lng, lat)
        if isfield(img.info.GPSInfo, {'GPSLatitude','GPSLongitude'})
          % Latitude (decimal degrees)
          lat = img.info.GPSInfo.GPSLatitude;
          lat = lat(1) + lat(2) / 60 + lat(3) / 3600;
          if deblank(strtrim(img.info.GPSInfo.GPSLatitudeRef)) == 'S'
            lat = -lat;
          end
          % Longitude (decimal degrees)
          lng = img.info.GPSInfo.GPSLongitude;
          lng = lng(1) + lng(2) / 60 + lng(3) / 3600;
          if deblank(strtrim(img.info.GPSInfo.GPSLongitudeRef)) == 'W'
            lng = -lng;
          end
          value.lnglat = [lng lat];
        end
        % GPS Altitude
        % NOTE: Unknown datum.
        if isfield(img.info.GPSInfo, 'GPSAltitude')
          value.altitude = img.info.GPSInfo.GPSAltitude; % (meters)
        end
        % GPS Satellites
        if isfield(img.info.GPSInfo, 'GPSSatellites')
          value.satellites = str2double(deblank(strtrim(img.info.GPSInfo.GPSSatellites)));
        end
      end
    end
    
    function value = get.size(img)
      value = [img.info.Width, img.info.Height];
    end
    
    function value = get.scale(img)
      if isempty(img.cam.imgsz)
        value = [];
      else
        value = Camera.getScaleFromSize(img.size, img.cam.imgsz);
      end
    end
    
    function set.scale(img, value)
      old_scale = img.scale;
      if isempty(old_scale)
        img.cam.imgsz = value * img.size;
      else
        img.cam = img.cam.resize(img.size);
        if value ~= 1
          img.cam = img.cam.resize(value);
        end
      end
      if ~isequal(old_scale, value)
        img.clear();
      end
    end

    % Image read
    
    function I = read(img, scale)
    % READ Read image data from file.
    % 
    %   I = img.read(scale = img.scale)
    % 
    % Inputs:
    %   scale - Scale factor for resizing the result
    
      if nargin < 2
        scale = img.scale;
      end
      if isempty(img.cached_I) || scale ~= img.scale
        I = imread(img.file);
        if ~isempty(scale) && scale ~= 1
          I = imresize(I, scale);
        end
        if scale == img.scale
          img.cached_I = I;
        end
      else
        I = img.cached_I;
      end
    end
    
    function clear(img)
    % CLEAR Clear object cache.
    % 
    %   img.clear()
    
      img.cached_I = [];
    end
  
  end % methods

  methods (Static)

  end % methods (Static)

  methods (Access = private)

  end % methods (Access = private)

end % classdef
