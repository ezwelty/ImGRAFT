function x = testf(varargin)

  p = inputParser;
  p.addOptional('x', []);
  p.addOptional('y', []);
  p.addOptional('z', []);
  p.parse(varargin);
  x = p.Results;

end
