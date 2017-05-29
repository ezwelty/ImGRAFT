function S = svg2struct(svg, scale)
  % SVG2STRUCT  Load all paths in a SVG XML file into a MATLAB structure.
  %
  %   S = svg2struct(svg, scale)
  %
  % Extracts the paths from an SVG file, grouping them according to layer
  % name in structure fields of the same name.
  %
  % NOTE: Only tested on SVG files generated by Adobe Illustrator CS6.
  %
  % Input:    svgfile    path to svg file
  %
  % Output:   S          structure with fields corresponding to svg layers
  %
  % See also svgpath2pts, img2mask, mask2pts.

  % Create XML DOM object
  if isempty(regexp(class(svg), 'org\.apache\.xerces\.dom\..*'))
    svg = xmlread(svg);
  end
  if ~strcmp(class(svg), 'org.apache.xerces.dom.DeepNodeListImpl') && svg.getElementsByTagName('svg').getLength > 0
    if nargin < 2
      viewbox = str2double(strsplit(char(svg.getElementsByTagName('svg').item(0).getAttribute('viewBox')), ' '));
      viewsz = [diff(viewbox([1 3])) diff(viewbox([2 4]))];
      images = svg.getElementsByTagName('image');
      if images.getLength > 0
        if images.getLength > 1
          warning('Multiple images found. Using first image to calculate scale factor.')
        end
        % TODO Use image transform tag instead? e.g. transform="matrix(0.24 0 0 0.24 -0.0015 0)"
        imgsz = [str2double(char(svg.getElementsByTagName('image').item(0).getAttribute('width'))), str2double(char(svg.getElementsByTagName('image').item(0).getAttribute('height')))];
        scale = imgsz ./ viewsz;
      else
        scale = [1 1];
      end
    else
      if length(scale) == 1
        scale = [scale scale];
      end
    end
    % Select Adobe Illustrator layers
    svg = svg.getElementsByTagName('svg');
  end

  % Initialize structure
  S = struct();

  % For each layer...
  for i = 0:(svg.getLength - 1)
    child = svg.item(i);
    tag = char(child.getNodeName);
    if ismember(tag, {'svg', 'switch'}) || (ismember(tag, {'g'}) && strcmp(child.getAttribute('id'), ''))
      S = svg2struct(child, scale);
      break
    end
    if ismember(tag, {'#comment', '#text', 'image', 'foreignObject'})
      continue
    end
    name = char(child.getAttribute('id'));
    name = regexprep(name, '_x([^_]+)_', '${char(hex2dec($1))}');
    if strcmp(name, '')
      name = tag;
    end
    fnames = fieldnames(S);
    if any(strcmp(name, fnames))
      nseries = sum(cell2mat(regexp(fnames, [name '_[0-9]+_'])));
      name = [name '_' num2str(nseries + 1) '_'];
    end
    switch tag
      case 'g'
        S.(name) = svg2struct(child, scale);
      case 'polyline'
        temp = char(child.getAttribute('points'));
        temp = textscan(temp, '%f', 'Delimiter', ' ,', 'MultipleDelimsAsOne', 1);
        temp = temp{1};
        data = reshape(temp, 2, length(temp) / 2)';
        S.(name) = bsxfun(@times, data, scale);
      case 'path'
        temp = char(child.getAttribute('d'));
        data = svgpath2pts(temp);
        S.(name) = bsxfun(@times, data, scale);
      case 'polygon'
        temp = char(child.getAttribute('points'));
        temp = textscan(temp, '%f%f', 'Delimiter', ',', 'MultipleDelimsAsOne', 1, 'EndOfLine', '\t');
        data = [temp{1} temp{2} ; temp{1}(1) temp{2}(1)];
        S.(name) = bsxfun(@times, data, scale);
      case 'line'
        x1 = str2num(child.getAttribute('x1'));
        y1 = str2num(child.getAttribute('y1'));
        x2 = str2num(child.getAttribute('x2'));
        y2 = str2num(child.getAttribute('y2'));
        data = [x1 y1; x2 y2];
        S.(name) = bsxfun(@times, data, scale);
      otherwise
        warning(['Unsupported tag encountered: ' tag])
    end
    if isstruct(S.(name)) && length(fieldnames(S.(name))) == 0
      S = rmfield(S, name);
    end
  end