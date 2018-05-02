function showRuler(ca)
axis manual
ruler = imline(ca);

% Get original position
pos = getPosition(ruler);

% Get updated position as the ruler is moved around
id = addNewPositionCallback(ruler,@(pos) title(sprintf('%s --x len: %s y len: %s ', mat2str(pos,3), num2str(pos(2,1)-pos(1,1)), num2str(pos(2,2)-pos(1,2)))));