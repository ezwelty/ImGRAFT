function sensorsz = sensorSize(mmakemodel)
  % SENSORSIZE  Get the CCD sensor width and height of a digital camera model.
  %
  %   sensorsz = sensorSize(makemodel)
  %
  % Returns the CCD sensor width and height in mm for the specified camera make
  % and model. Data is from Digital Photography Review (www.dpreview.com).
  % See also https://www.dpreview.com/articles/8095816568/sensorsizes.
  %
  % Inputs:
  %   makemodel - camera make and model [make ' ' model]
  %
  % Outputs:
  %   sensorsz - sensor size in mm [width height]

  % Check number of arguments
  if (nargin < 1)
      error('Specify make & model of camera.')
  end

  % Check type of arguments
  if (~ischar(makemodel))
      error('Camera make & model must be a character string: [Make Model]');
  end

  % Data
  % Camera names: [deblank(strtrim(Make)) ' ' deblank(strtrim(Model))]
  % Sensor sizes (mm)
  data = {
    'NIKON CORPORATION NIKON D2X', [23.7 15.7]; % https://www.dpreview.com/reviews/nikond2x/2
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
