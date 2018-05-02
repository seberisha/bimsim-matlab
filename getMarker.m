function m = getMarker(i)
markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
m = markers{mod(i,numel(markers))+1};