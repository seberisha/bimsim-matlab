function x = changeRange(y, NewMin, NewMax)
OldMax = max(y(:));
OldMin = min(y(:));
if isscalar(y)
    OldMin = 0;
end
OldRange = (OldMax - OldMin);  

% if (OldRange == 0)
%     OldRange = OldMax;
% end


NewRange = (NewMax - NewMin); 
x = (((y - OldMin) .* NewRange) ./ OldRange) + NewMin;

% OldRange = (OldMax - OldMin);
% if (OldRange == 0)
%     NewValue = NewMin;
% else
%     NewRange = (NewMax - NewMin)  ;
%     x = (((y - OldMin) * NewRange) / OldRange) + NewMin;
% end