function hh = scrollView(vol, dim, limits)


app=[];

app.size = size(vol);

if min(app.size) == 1
	fprintf('Wrong size: ');
	disp(app.size);
	return
end

if ~isreal(vol)
	error('Only real values may displayed!')
end

if nargin > 2
	if limits == 0
		app.limits = [min(vol(:)), max(vol(:))];
		if app.limits(1) == app.limits(2)
			app.limits(1) = app.limits(1) -eps;
			app.limits(2) = app.limits(1) +eps;
		end
	else
		app.limits = limits;
	end
else
	app.limits = [];
end

if nargin < 2
	app.dim = 3;
else
	app.dim = dim;
end

maxval = size(vol, app.dim);

if nargout > 0
	hh = figure;
else
	figure
end
%imagesc( reshape( vol(1000:1013:1013*1013*601), 1013, [])')
    colormap jet(1024)
  %  colormap(gray);

app.axis = gca;
app.figure = gcf;

% Generate constants for use in uicontrol initialization
pos=get(app.axis,'position');
Newpos=[pos(1) pos(2)-0.1 pos(3) 0.05];


% Creating Uicontrol
h=uicontrol('style','slider',...
    'units','normalized','position',Newpos,...
    'callback',@sfScroll,'min',1, 'max', maxval, 'value', floor(maxval/2), 'sliderstep', [1/maxval 10/maxval] );

%imagesc(zeros(50));
sfScroll( h, [] );

%%%%%%%%%%%%%%%%%%%%%% subfunctions %%%%%%%%%%%%%
% 
function sfScroll(varargin)

handle = varargin{1};
%scrollVal = floor(get(gcbo,'value'));
scrollVal = floor(get(handle,'value'));

	switch app.dim
		case 1
			sfImagescWrapper( vol(scrollVal, :, :), app.limits )
		case 2
			sfImagescWrapper( vol(:, scrollVal, :), app.limits )
		case 3
			sfImagescWrapper( vol(:, :, scrollVal), app.limits )
	end

	title(scrollVal);
	
end


function sfImagescWrapper(image, limits)
	
	if ( length(limits) >= 2 && limits(1) ~= limits(2))
		imagesc( squeeze(image), limits)

	else
		imagesc( squeeze(image))


	end

	axis image
	colorbar
end



end
