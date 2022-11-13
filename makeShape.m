function xy = makeShape(shape, L, location, w, h, theta)

if nargin < 6,
    theta = 0;
end;


xC = w/2 + location(1)*w; % x position of center of shape, in pixels
yC = h/2 + location(2)*h; % y position of center of shape, in pixels

switch shape,
    case 'square',
        L = L(1);
        xy = [- L/2, - L/2; ...
              + L/2, - L/2; ...
              + L/2, + L/2; ...
              - L/2, + L/2];
    case 'rectangle',
        xy = [- L(1)/2, - L(2)/2; ...
              + L(1)/2, - L(2)/2; ...
              + L(1)/2, + L(2)/2; ...
              - L(1)/2, + L(2)/2];
    otherwise,
        fprintf('makeShape.m: don"t know how to make shape %s, so returning a square\n', shape);
        L = L(1);
        xy = [- L/2, - L/2; ...
              + L/2, - L/2; ...
              + L/2, + L/2; ...
              - L/2, + L/2];
              
end;

if abs(theta) > eps,
    % rotations about the center of the shape
    A = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    xy = (A*xy')';
end;

xy(:,1) = xC + xy(:,1);
xy(:,2) = yC + xy(:,2);