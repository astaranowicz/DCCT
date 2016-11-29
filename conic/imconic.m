function [kind,a,b,R,t] = imconic(p, figax, ellipseColor, w,varargin)
% Compute parameters of or plot conic section given in implicit form.
% The conic section is assumed to be given as
% a*x^2 + b*x*y + c*y^2 + d*x + f*y + g = 0
% where [a,b,c,d,f,g] are stacked into the parameter vector p.
% If given, only those points of the conic are plotted that fall within
% the specified window.
%
% Input arguments:
% p:
%    the parameter vector p = [a b c d f g]
% w:
%    the window w = [x_min x_max y_min y_max], or 0 to suppress plotting
%    even for call with single output argument
% figax:
%    axes to use for plotting the conic section
%
% Added by AS:
% ellipseColor - change the color of the ellipse
% intolerance - tolerance on the functions such as a-c <tol
%
%
% Output arguments:
% h:
%    handle of the graphics object created
% -- OR --
% kind:
%    the type of the conic section: circle, ellipse, parabola, hyperbola
% a, b:
%    explicit parameters (semi-axes) of conic section
% R:
%    rotation from standard alignment
% t:
%    offset from origin
%
% Examples:
% h = imconic([4 0 4 4 -36 -16])
%    plots the conic section in the current axes, returning its handle h
% h = imconic([4 0 4 4 -36 -16], [], a)
%    plots the conic section in axes a, returning its handle h
% [kind,a,b,R,t] = imconic([1 0 0 0 1 0])
%    computes parameters for the conic section
% kind = imconic([4 0 4 4 -36 -16], 0)
%    identifies but does not plot the conic section

% Copyright 2010 Levente Hunyadi
%%
% global ellipseColor2;
% if ~exist('ellipseColor'),
%     ellipseColor2 = 'b';
% else
%     ellipseColor2 = ellipseColor;
% end
%%
validateattributes(p, {'numeric'}, {'nonempty','real','vector'});
assert(length(p) == 6, 'imconic:DimensionMismatch', ...
    'A parameter vector of length 6 expected for a to g in a*x^2 + b*x*y + c*y^2 + d*x + f*y + g = 0.');
    tolerance2 = 1e-12;
%figure that will display the plot
if nargin > 1 && ~isempty(figax)
    validateattributes(figax, {'numeric'}, {'scalar'});
    assert(ishandle(figax), 'imconic:InvalidArgumentValue', ...
        'The figure axes parameter is expected to be a valid graphics handle.');
else
    figax = [];
end
% the color of the drawn conic
if nargin > 2 && ~isempty(ellipseColor)
    ellipseColor2 = ellipseColor;
else
    ellipseColor2 = 'b';
end
% window size
if nargin > 3 && ~isempty(w)
    if isscalar(w)
        assert(w == 0, 'imconic:InvalidArgumentValue', ...
            'A parameter vector of length 4 expected as window, or use 0 to suppress plotting.');
    else
        validateattributes(w, {'numeric'}, {'nonempty','real','vector'});
        assert(numel(w) == 4, 'imconic:DimensionMismatch', ...
            'A parameter vector of length 4 expected for x_min, x_max, y_min and y_max.');
    end
else
    w = [];
end


tol = tolerance2;
tol_discr_delta = 1e-12;
% tol_abc = 1e-1;
% conic section discriminants
J = imconicdiscr(p);
Delta = imconicdelta(p);
t = imconictranslation(p); %Note: before the normalization
p = p/abs(p(6));

if nargout > 1 || isempty(w) || ~isscalar(w) || w ~= 0  % 2 or more output arguments, or plot conic section
    % choose kind of conic section
%     if abs(J) > tol_discr_delta,  % equation represents an ellipse
%         a = p(1);
%         b = 0.5 * p(2);
%         c = p(3);
        %assert(Delta ~= 0, 'imconic:InvalidArgumentValue', ...
        %    'Parameters define a degenerate conic: imaginary intersecting lines.');
        %assert(Delta / (a+c) < 0, 'imconic:InvalidArgumentValue', ...
        %    'Parameters define an imaginary ellipse.');
        %         if a == c && b == 0  % equation represents a circle
%         if (abs(a-c) < tol) && (abs(b) < tol),
%             kind = 'circle';
%             if sign(-p(6) + p(1)*t(1)^2 + p(2)*t(1)*t(2) + p(3)*t(2)^2) == -1
%                 a = sqrt(-p(6) - p(1)*t(1)^2 + p(2)*t(1)*t(2) - p(3)*t(2)^2);
%             else
%                 a = realsqrt(t(1)^2 + t(2)^2 - ((-2*t(1))/p(4)));
%                 a = sqrt(-p(6) + p(1)*t(1)^2 + p(2)*t(1)*t(2) + p(3)*t(2)^2);
%             end
            %a = imcircle(p);
%             b = a;  % circle radius
            %display('circle');
%         else
%             p_t = p/abs(p(6));
            kind = 'ellipse';
            %display('ellipse');
            [a,b] = imellipse(p);
%         end
        %elseif J < -1e-12  % equation represents a hyperbola
        %    assert(Delta ~= 0, 'imconic:InvalidArgumentValue', ...
        %        'Parameters define a degenerate conic: real intersecting lines.');
        %    kind = 'hyperbola';
        %    [a,b] = imhyperbola(p);
%     else  % equation represents a parabola
%         if abs(Delta) < tol_discr_delta,
%             a = p(1);
%             c = p(3);
%             d = 0.5 * p(4);
%             f = 0.5 * p(5);
%             g = p(6);
%             K = det([a d ; d g]) + det([c f ; f g]);
%             if K > 0
%                 error('imconic:InvalidArgumentValue', ...
%                     'Parameters define a degenerate conic: imaginary parallel lines.');
%             elseif K < 0
%                 error('imconic:InvalidArgumentValue', ...
%                     'Parameters define a degenerate conic: real parallel lines.');
%             else
%                 error('imconic:InvalidArgumentValue', ...
%                     'Parameters define a degenerate conic: coincident lines.');
%             end
%         end
%         kind = 'parabola';
%         a = imparabola(p);
%         b = NaN;  % meaningless for parabola
%     end
    
    %% Rotation calculation
%     p = p/abs(p(6));
%     if strcmp(kind, 'circle') == 1, % if so, then it is a circle
%         R = [1 0; 0 1];
%     else %Is an ellipse
        C_2x2 = [  p(1)    p(2)/2;
            p(2)/2   p(3) ];
        % U is in the form [-s c; c s]
        [Ur,Dr,Vr] = svd(C_2x2);
        % Since 0<=alpha<=pi ==> sin>=0 ==> U(2,2)>0
        Ur = Ur*sign(Ur(2,2));
        Vr = Vr*sign(Ur(2,2));
        % Rotation
        if Dr(1,1)~=Dr(2,2),
            R = [Ur(:,2) Ur(:,1)];
        else
            R = Ur;
        end
%     end
    
    %R = imconicrotation(imconictranslate(p, -t));
    
    if nargout < 2  % h = imconic(p)
        zx = imconicisect(p,w,t,R);
        
        switch kind
            case 'circle'
                [z,grp] = imcirclepoints(a,t,zx,w);
            case 'ellipse'
                [z,grp] = imellipsepoints(a,b,R,t,zx,w);
            case 'parabola'
                [z,grp] = imparabolapoints(a,R,t,zx,w);
            case 'hyperbola'
                [z,grp] = imhyperbolapoints(a,b,R,t,zx,w);
        end
        kind = implotconic(figax,z,grp,R,t,ellipseColor2,varargin{:});
    end
else  % kind = imconic(p,0)
    % choose kind of conic section
    if J > tol  % equation represents an ellipse
        a = p(1);
        b = 0.5 * p(2);
        c = p(3);
        if Delta == 0
            kind = 'imaginary_intersecting_lines';
        elseif Delta / (a+c) >= 0
            kind = 'imaginary_ellipse';
        elseif a == c && b == 0
            kind = 'circle';
        else
            kind = 'ellipse';
        end
    elseif J < -tol  % equation represents a hyperbola
        if Delta == 0
            kind = 'real_intersecting_lines';
        else
            kind = 'hyperbola';
        end
    else  % equation represents a parabola
        if Delta == 0
            kind = 'parabola';
        else
            a = p(1);
            c = p(3);
            d = 0.5 * p(4);
            f = 0.5 * p(5);
            g = p(6);
            K = det([a d ; d g]) + det([c f ; f g]);
            if K > 0
                kind = 'imaginary_parallel_lines';
            elseif K < 0
                kind = 'real_parallel_lines';
            else
                kind = 'coincident_lines';
            end
        end
    end
end

%%
function Delta = imconicdelta(p)

A_xx = p(1);
A_xy = 0.5*p(2);
A_yy = p(3);
B_x  = 0.5*p(4);
B_y  = 0.5*p(5);
C    = p(6);
Delta = det([ A_xx A_xy B_x ;
    A_xy A_yy B_y ;
    B_x  B_y  C ]);
%%
function h = implotconic(figax, z, grp, R, t,ellipseColor2,varargin)
% Plot conic section.
% Points are passed in groups, each of whose points are connected with lines,
% but the same color and marker is used for plotting the conic section.
%
% Input arguments:
% figax:
%   the graphics handle of the axes to use for plotting
% z:
%    a 2-by-n matrix of points to plot
% grp:
%    a vector of group indices for each point in z
% R:
%    the rotation matrix
% t:
%    the translation vector
%
% Output arguments:
% h:
%    handle of the graphics object created
% global ellipseColor2;
if isempty(z)
    return;  % nothing to plot
end
if isempty(figax)
    figax = gca;
end

z = transrot(z, t, R);
x = z(1,:);
y = z(2,:);

fig = ancestor(figax, 'figure');
fignextplot = get(fig, 'NextPlot');   % save setting for figure
axnextplot = get(figax, 'NextPlot');  % save setting for axes
set([fig,figax], 'NextPlot', 'add');  % re-use same figure

% get color and line style of first segment
grpcount = max(grp);
% if grpcount > 1
%     h = hggroup('Parent', figax);
%
%     lineseries = plot(figax, x(grp==1), y(grp==1), 'Parent', h, varargin{:});
%     color = get(lineseries, 'Color');
%     linestyle = get(lineseries, 'LineStyle');
%     marker = get(lineseries, 'Marker');
%
%     for k = 2 : grpcount
%         plot(figax, x(k==grp), y(k==grp), ...
%             'Parent', h, ...
%             'Color', color, ...  % re-use color and line style for other line segments
%             'LineStyle', linestyle, ...
%             'Marker', marker, ...
%             varargin{:});
%     end
% else
h = plot(gca, x(grp==1), y(grp==1), varargin{:}, 'Color', ellipseColor2);
% end

set(fig, 'NextPlot', fignextplot);   % reset figure re-use
set(figax, 'NextPlot', axnextplot);  % reset axes re-use
%%
function r = imcircle(p)
% Compute explicit parameters for a circle from its implicit equation.
% This function uses the implicit circle equation
% a*x^2 + a*y^2 + 2*d*x + 2*f*y + g = 0
% where p = [a 0 c 2*d 2*f g].
%
% Output arguments:
% t:
%    translation vector
% r:
%    circle radius

validateattributes(p, {'numeric'}, {'nonempty','real','vector'});
% assert(length(p) == 6, 'imconic:DimensionMismatch', ...
%     'A parameter vector of length 6 expected for a to g in a*x^2 + b*x*y + c*y^2 + d*x + f*y + g = 0.');
% assert(p(2) == 0 && p(1) == p(3), 'imconic:InvalidArgumentValue', ...
%     'Invalid parameter vector for a circle, expected b == 0 and a == c.');

% extract parameters from vector form
a = p(1);
d = p(4);
e = p(5);
f = p(6);

if (f) == -1
    r = realsqrt(((d^2+e^2)/(4*a^2))+(1/a));
else
    % use implicit equation to compute circle radius
    r = realsqrt(((d^2+e^2)/(4*a^2))-(1/a));
end



%%
function [z,grp] = imcirclepoints(r, t, zx, w)
% Plot a circle given its explicit parameters.
%
% Input arguments:
% t:
%    translation vector
% r:
%    circle radius

[alpha,grp] = implotangles(atan2(zx(2,:), zx(1,:)), -pi, [-r;0], w, t, eye(2,2));
x = r .* cos(alpha);
y = r .* sin(alpha);
z = [ x ; y ];
%%
function [semi_a,semi_b] = imellipse(p)
% Compute explicit parameters for an ellipse from its implicit equation.
% The function uses the implicit ellipse equation:
% a*x^2 + 2*b*x*y + c*y^2 + 2*d*x + 2*f*y + g = 0
% where p = [a 2*b c 2*d 2*f g].

validateattributes(p, {'numeric'}, {'nonempty','real','vector'});
% assert(length(p) == 6, 'imconic:DimensionMismatch', ...
%     'A parameter vector of length 6 expected for a to g in a*x^2 + b*x*y + c*y^2 + d*x + f*y + g = 0.');
% assert(imconicdiscr(p) > 0, 'imconic:InvalidArgumentValue', ...
%     'Invalid parameter vector for an ellipse, expected b^2 - 4*a*c < 0.');
% assert(imconicdelta(p) ~= 0, 'imconic:InvalidArgumentValue', ...
%     'Parameters define a degenerate ellipse: imaginary intersecting lines.');

% extract parameters from vector form
a = p(1);
b = 0.5 * p(2);
c = p(3);
d = 0.5 * p(4);
f = 0.5 * p(5);
g = p(6);

% use implicit equation to compute ellipse semi-axes length
q = 2*(a*f^2+c*d^2+g*b^2-2*b*d*f-a*c*g)/(b^2-a*c);
r = realsqrt((a-c)^2+4*b^2);

if sign(q/(r-(a+c))) == -1
    keyboard
elseif sign(q/(-r-(a+c))) == -1
    keyboard
end
semi_a = realsqrt(q/(r-(a+c)));   % major axis
semi_b = realsqrt(q/(-r-(a+c)));  % minor axis

%%
function [z,grp] = imellipsepoints(a, b, R, t, zx, w)
% Plot an ellipse given its explicit parameters.
%
% Input arguments:
% R:
%    rotation matrix to transform the ellipse in the standard coordinate
%    system into the target system
% t:
%    translation vector to transform the ellipse in the standard coordinate
%    system into the target system
% zx:
%    coordinates of points in the standard system where the conic and the
%    window intersect
% w:
%    window boundaries in the standard coordinate system

[alpha,grp] = implotangles(atan2(a/b .* zx(2,:), zx(1,:)), -pi, [-a;0], w, t, R);
x = a.*cos(alpha);
y = b.*sin(alpha);
z = [ x ; y ];
%%
function [a,b] = imhyperbola(p)
% Compute explicit parameters for a hyperbola from its implicit equation.
% The function uses the implicit hyperbola equation:
% A_xx*x^2 + 2*A_xy*x*y + A_yy*y^2 + 2*B_x*x + 2*B_y*y + C = 0
% where p = [A_xx 2*A_xy A_yy 2*B_x 2*B_y C].

validateattributes(p, {'numeric'}, {'nonempty','real','vector'});
assert(length(p) == 6, 'imconic:DimensionMismatch', ...
    'A parameter vector of length 6 expected for a to g in a*x^2 + b*x*y + c*y^2 + d*x + f*y + g = 0.');
D = imconicdiscr(p);
assert(D < 0, 'imconic:InvalidArgumentValue', ...
    'Invalid parameter vector for a hyperbola, expected b^2 - 4*a*c > 0.');
Delta = imconicdelta(p);
assert(Delta ~= 0, 'imconic:InvalidArgumentValue', ...
    'Parameters define a degenerate hyperbola: real intersecting lines.');

A_xx = p(1);
A_yy = p(3);

% major and minor semi-axes
lambda = roots([1, -(A_xx+A_yy), D]);  % lambda^2 - (A_xx + A_yy)*lambda + D = 0
a = realsqrt(abs(Delta/(lambda(1)*D)));
b = realsqrt(abs(Delta/(lambda(2)*D)));
%%
function [z,grp] = imhyperbolapoints(a, b, R, t, zx, w)
% Plot a hyperbola given its explicit parameters.

r = a/b .* (zx(2,:) ./ zx(1,:));
r(r > 1 | r < -1) = NaN;
alphax = atanh(r);

% left branch
xlx = -a.*cosh(alphax);
ylx = -b.*sinh(alphax);
f = abs(zx(1,:) - xlx) < 1e-12 & abs(zx(2,:) - ylx) < 1e-12;  % find points that are on left branch
[alpha,grpl] = implotangles(alphax(f), 0, [-a;0], w, t, R);
xl = -a.*cosh(alpha);
yl = -b.*sinh(alpha);

% right branch
xrx = a.*cosh(alphax);
yrx = b.*sinh(alphax);
f = abs(zx(1,:) - xrx) < 1e-12 & abs(zx(2,:) - yrx) < 1e-12;  % find points that are on right branch
[alpha,grpr] = implotangles(alphax(f), 0, [a;0], w, t, R);
xr = a.*cosh(alpha);
yr = b.*sinh(alpha);

x = [ xl, xr ];
y = [ yl, yr ];
z = [ x ; y ];

if ~isempty(grpl) && ~isempty(grpr)
    grp = [ grpl, max(grpl)+grpr ];
elseif ~isempty(grpl)
    grp = grpl;
elseif ~isempty(grpr)
    grp = grpr;
else
    grp = [];
end
%%

function a = imparabola(p)

validateattributes(p, {'numeric'}, {'nonempty','real','vector'});
% assert(length(p) == 6, 'imconic:DimensionMismatch', ...
%     'A parameter vector of length 6 expected for a to g in a*x^2 + b*x*y + c*y^2 + d*x + f*y + g = 0.');
% D = imconicdiscr(p);
% assert(abs(D) < 1e-3, 'imconic:InvalidArgumentValue', ...
%     'Invalid parameter vector for a parabola, expected b^2 = 4*a*c.');

[R,theta] = imconicrotation(p);
p = imconicrotate(p, -theta);

if abs(p(4)) > 1e-12
    % x = a*y^2 + b*y + c, p = [0 0 a 1 b c]
    p = p ./ p(4);
    a = p(3); %b = p(5); c = p(6);
    %x0 = (4*a*c - b^2) / (4*a);
    %y0 = -b / (2*a);
elseif abs(p(5)) > 1e-12
    % y = a*x^2 + b*x + c, p = [a 0 0 b 1 c]
    p = p ./ p(5);
    a = p(1); %b = p(4); c = p(6);
    %x0 = -b / (2*a);
    %y0 = (4*a*c - b^2) / (4*a);
end
%%
function [z,grp] = imparabolapoints(a, R, t0, zx, w)

[t,grp] = implotangles(2*zx(1,:) ./ zx(2,:), 0, [0;0], w, t0, R);
x = a.*t.^2;
y = 2*a.*t;
z = [ x ; y ];
%%
function tf = isinrect(z, w, t, R)
% Points that are inside a window.
%
% Input arguments:
% z:
%    a 2-by-n matrix of points
% w:
%    window boundaries as [x_min x_max y_min y_max]
% t:
%    the translation vector applied to the points
% R:
%    the rotation matrix applied to the points

if nargin > 3
    z = transrot(z, t, R);
elseif nargin > 2
    z = transrot(z, t);
end

x = z(1,:); y = z(2,:);
x_min = w(1); x_max = w(2);
y_min = w(3); y_max = w(4);

tf = x >= x_min & x_max >= x & y >= y_min & y_max >= y;
%%
function [alpha,grp] = implotangles(alphax, beta, z, w, t, R)
% Input arguments:
% alphax:
%    polar coordinate angles of intersections with bounds window [rad]
% z:
%    a point that is part of the conic
% beta:
%    polar coordinate angle of point z in the standard coordinate system
%
% Output arguments:
% alpha:
%    polar coordinate angles to plot [rad]
% bnd:
%    angle bounds [rad]

alpha = linspace(-pi, pi, 200);
if ~isempty(w)
    alpha = sort([ alphax-1e-12 alphax+1e-12 alpha ]);
    alphax = sort(alphax);
    bnd = reshape(alphax, 2, numel(alphax)/2);
    if isinrect(z, w, t, R) == isempty(intfilter(beta, bnd))  % point of conic is inside window but discarded by filter or vice versa
        alphax = [-pi alphax pi];  % swap inside and outside
        bnd = reshape(alphax, 2, numel(alphax)/2);
    end
    [alpha,grp] = intfilter(alpha, bnd);
else
    grp = ones(1,numel(alpha));
end
%%