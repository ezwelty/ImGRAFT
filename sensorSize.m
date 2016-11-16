function sensorsz = sensorSize(varargin)
  % SENSORSIZE  Get the CCD sensor width and height of a digital camera model.
  %
  %   sensorsz = sensorSize(makemodel)
  %   sensorsz = sensorSize(make, model)
  %
  % Returns the CCD sensor width and height in mm for the specified camera make
  % and model. Data is from Digital Photography Review (www.dpreview.com).
  % See also https://www.dpreview.com/articles/8095816568/sensorsizes.
  %
  % Inputs:
  %   makemodel - camera make and model [make ' ' model]
  %   make      - camera make
  %   model     - camera model
  %
  % Outputs:
  %   sensorsz - sensor size in mm [width height]

  % Check input arguments
  if (nargin < 1)
      error('Specify make & model of camera.')
  end
  if (nargin > 1)
    makemodel = [deblank(strtrim(varargin{1})) ' ' deblank(strtrim(varargin{2}))];
  end

  % Sensor sizes (mm)
  data = {
    'NIKON CORPORATION NIKON D2X', [23.7 15.7]; % https://www.dpreview.com/reviews/nikond2x/2
    'NIKON CORPORATION NIKON D200', [23.6 15.8]; % https://www.dpreview.com/reviews/nikond200/2
  };

  % Check for match
  match = find(strcmp(makemodel, data(:, 1)));

  % If match found, return sensor size
  if (~isempty(match))
      sensorsz = data{match, 2};
  else
      warning(['No sensor size found for camera "' makemodel '". Returning [].']);
      sensorsz = [];
  end
