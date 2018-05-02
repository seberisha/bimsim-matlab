function imageExample()
    %# RGB image
    img = imread('peppers.png');
    sz = size(img);

    %# show image
    hFig = figure();
    hAx = axes();
    image([1 sz(2)]+100, [1 sz(1)]+200, img)    %# shifted XData/YData

    %# hook-up mouse button-down event
    set(hFig, 'WindowButtonDownFcn',@mouseDown)

    function mouseDown(o,e)
        %# get current point
        p = get(hAx,'CurrentPoint');
        p = p(1,1:2);

        %# convert axes coordinates to image pixel coordinates
        %# I am also rounding to integers
        x = round( axes2pix(sz(2), [1 sz(2)], p(1)) );
        y = round( axes2pix(sz(1), [1 sz(1)], p(2)) );

        %# show (x,y) pixel in title
        title( sprintf('image pixel = (%d,%d)',x,y) )
    end
end
