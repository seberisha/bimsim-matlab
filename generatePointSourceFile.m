function ps = generatePointSourceFile(thresh, rows, cols)
% %lambda = 2.543;
% % step_r = fov/rows;
% % step_c = fov/cols;
% % space_r = round(lambda/step_r);
% % space_c = round(lambda/step_c);
% % 
% % if (space_r < 1)
% %     space_r = 1;
% % end
% % 
% % if(space_c <1)
% %     space_c=1;
% % end
% % ps = zeros(rows,cols);
% % ps(1:space_r:end, 1:space_c:end) = 255;

% ps = zeros(rows,cols);
% ps(1:floor(fov/lambda):end,1:floor(fov/lambda):end) = 255;

[rr, cc] = meshgrid(1:rows,1:cols);
ps = sqrt((rr-rows/2).^2+(cc-cols/2).^2)<=thresh;
ps = ps.*255;
figure, imshow(ps)
%imwrite(ps, sprintf('%iPointSources.png',numel(find(ps==255))))