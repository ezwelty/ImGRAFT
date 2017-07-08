function [] = flip_images(img1, img2, cols, rows, times)
  if nargin < 3 || isempty(cols)
    cols = 1:size(img1, 2);
    rows = 1:size(img1, 1);
  end
  if nargin == 3
    rows = cols;
  end
  if nargin < 5
    times = [];
  elseif times == 'auto'
    times = Inf;
  end
  figure('units', 'normalized', 'outerposition', [0 0 1 1]);
  i = 0;
  while true
    if i == 0
      h = image(img1(rows, cols, :));
      title('1');
    elseif mod(i, 2) == 0
      set(h, 'cdata', img1(rows, cols, :));
      title('1');
    else
      set(h, 'cdata', img2(rows, cols, :));
      title('2');
    end
    if isempty(times)
      waitforbuttonpress;
    elseif times == i
      break;
    else
      pause(0.25);
    end
    i = i + 1;
  end
end
