function P = cpuLegendre(order,x)
%P(:,:,1) = ones(size(x));
% if(order==0)
%     return;
% end
% if isscalar(x)
%     
%     P(2) = x;
%     for j=2:order
%         %matlab definition
%         P(j+1) = ((2*j-1)/j).*x.*(P(j)) - ((j-1)/j).*(P(j-1));
%         %other definition
%         %P(:,:,j+1) = ((2*j+1)/(j+1)).*x.*squeeze(P(:,:,j))- ((j)/(j+1)).*squeeze(P(:,:,j-1));
%     end
% else
%     
%     P(:,:,2) = x;
%     for j=2:order
%         %matlab definition
%         P(:,:,j+1) = ((2*j-1)/j).*x.*(P(:,:,j)) - ((j-1)/j).*(P(:,:,j-1));
%         %other definition
%         %P(:,:,j+1) = ((2*j+1)/(j+1)).*x.*squeeze(P(:,:,j))- ((j)/(j+1)).*squeeze(P(:,:,j-1));
%     end
% end



if isscalar(x)
    P = zeros(order+1,1);
    P(1) = 1;
    if(order==0)
        return;
    end
    P(2) = x;
    for j=2:order
        %matlab definition
        P(j+1) = ((2*j-1)/j).*x.*(P(j)) - ((j-1)/j).*(P(j-1));
        %other definition
        %P(:,:,j+1) = ((2*j+1)/(j+1)).*x.*squeeze(P(:,:,j))- ((j)/(j+1)).*squeeze(P(:,:,j-1));
        
    end
elseif isvector(x)
    P = zeros(numel(x),order+1);
    P(:,1) = 1;
    if(order==0)
        return;
    end
    P(:,2) = x;
    for j=2:order
        %matlab definition
        P(:,j+1) = ((2*j-1)/j).*x.*(P(:,j)) - ((j-1)/j).*(P(:,j-1));
        %other definition
        %P(:,:,j+1) = ((2*j+1)/(j+1)).*x.*squeeze(P(:,:,j))- ((j)/(j+1)).*squeeze(P(:,:,j-1));
        
    end
else
    P = zeros([size(x) order+1]);
    P(:,:,1) = 1;
    if(order==0)
        return;
    end
    P(:,:,2) = x;
    for j=2:order
        %matlab definition
        P(:,:,j+1) = ((2*j-1)/j).*x.*(P(:,:,j)) - ((j-1)/j).*(P(:,:,j-1));
        %other definition
        %P(:,:,j+1) = ((2*j+1)/(j+1)).*x.*squeeze(P(:,:,j))- ((j)/(j+1)).*squeeze(P(:,:,j-1));
        
    end
end
