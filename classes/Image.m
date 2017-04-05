classdef Image
  % Image Photographic image data structure
  %
  % Image Properties:
  % imgsz    - Size of image in pixels [nx|ncols|width, ny|nrows|height]
  % f        - Focal length in pixels [fx, fy]
  % c        - Camera center coordinates in pixels [cx, cy]
  % k        - Radial distortion coefficients [k1, ..., k6]
  % p        - Tangential distortion coefficients [p1, p2]
  % xyz      - Camera position in world coordinates [x, y, z]
  % viewdir  - Camera view direction in degrees [yaw, pitch, roll]
  %            yaw: clockwise rotation about z-axis (0 = look north)
  %            pitch: rotation from horizon (+ look up, - look down)
  %            roll: rotation about optical axis (+ down right, - down left, from behind)
  % sensorsz - Size of camera sensor in mm [width, height] (optional)
  %
  % Image Properties (dependent):
  % fullmodel - Vector containing all 20 camera parameters
  %             [xyz(1:3), imgsz(1:2), viewdir(1:3), f(1:2), c(1:2), k(1:6), p(1:2)]
  % fmm       - Focal length in mm [fx, fy] (not set unless sensorsz is defined)
  % R         - Rotation matrix corresponding to camera view direction (read-only)
  % K         - Camera matrix [fx 0 cx; 0 fy cy; 0 0 1] (read-only)
  %
  % Image Methods:
  % Image      - Construct a new Image

  properties
    info
    cam
  end

  properties (Dependent)
    file
    datestr
    datenum
    shutter
    aperture
    iso
    ev
    gps
  end

  methods

    % Image creation

    function img = Image(filename, cam)
      if nargin > 1 && class(cam) ~= 'camera'
        error('Not an object of class camera.');
      end
      % Metadata
      img.info = imfinfo(filename);
      % Image size in pixels [nx, ny]
      if nargin > 1
        if isempty(cam.imgsz)
          imgcam = cam;
          imgcam.imgsz = [img.info.Width img.info.Height];
        else
          imgcam = cam.resize([img.info.Width img.info.Height]);
        end
      else
        imgcam = camera('imgsz', [img.info.Width img.info.Height]);
      end
      % Sensor size [mm width, mm height]
      if isempty(imgcam.sensorsz)
        if isfield(img.info, {'Make', 'Model'})
          imgcam.sensorsz = sensorSize(img.info.Make, img.info.Model);
        end
      end
      % Focal length in mm (if f not already set)
      if ~isempty(imgcam.sensorsz) && isempty(imgcam.f)
        if isfield(img.info.DigitalCamera, 'FocalLength')
          imgcam.fmm = img.info.DigitalCamera.FocalLength;
        elseif isfield(img.info.DigitalCamera, 'FocalLengthIn35mmFilm')
          % FIXME: Convert to true focal length using sensor size?
          imgcam.fmm = img.info.DigitalCamera.FocalLengthIn35mmFilm;
          warning('True focal length not found, using 35mm equivalent.');
        end
      end
      % Camera
      img.cam = imgcam;
    end

    % Image properties

    function value = get.file(img)
      value = img.info.Filename;
    end

    function value = get.datestr(img)
      if isfield(img.info.DigitalCamera, 'DateTimeOriginal')
        value = deblank(img.info.DigitalCamera.DateTimeOriginal);
        if isfield(img.info.DigitalCamera, 'SubsecTimeOriginal')
          value = [value '.' deblank(strtrim(img.info.DigitalCamera.SubsecTimeOriginal))];
        elseif isfield(img.info.DigitalCamera, 'SubsecTime')
          value = [value '.' deblank(strtrim(img.info.DigitalCamera.SubsecTime))];
        end
      else
        value = [];
      end
    end

    function value = get.datenum(img)
      value = datenum(img.datestr, 'yyyy:mm:dd HH:MM:SS.FFF');
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

    % Image read

    function I = read(img)
      I = imread(img.file);
    end

  end % methods

  methods (Static)

  end % methods (Static)

  methods (Access = private)

  end % methods (Access = private)

end % classdef
