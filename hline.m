function ax = hline(x,varargin)
%% hline(x)
% Plots horizontal line on current axes at X=x. (x can be a vector)
% Default is a black dashed line, but user can specify formatting

% if first input is axes
if isnumeric(x),
    ax = gca;
elseif strcmp(get(x, 'type'), 'axes'),
    ax = x;
    x = varargin{1};
    varargin = varargin(2:end);
else
    error('invalid input')
end

x = x(:)';
X = [repmat(x,2,1);nan(2,length(x))];
X = X(:);

yy = get(ax,'XLim')';
Y = [repmat(yy,1,length(x));nan(2,length(x))];
Y = Y(:);

prop = get(ax,'NextPlot');
set(ax,'NextPlot','add')
if nargin < 2,
    plot(ax,Y,X,'--k');
else
    plot(ax,Y,X,varargin{:});
end
set(ax,'NextPlot',prop)