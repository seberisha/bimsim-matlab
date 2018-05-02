function [wav, absImages, spec] = loadFTIR(hdrPath, absImgPath, add, x_pos, y_pos)
%load data from Cary FTIR

addpath(genpath('/home/sberisha/source/stimlib/stim/matlab'))
eh = rtsEnviLoadHeader(hdrPath);
wav = eh.wavelength;
absImages = multibandread(absImgPath,[eh.lines,eh.samples,eh.bands],'float',eh.header_offset, eh.interleave,'ieee-le');

%absImages(absImages<0)=0;
spec = zeros(size(wav));
numWav = numel(wav);
r_idx = round(eh.lines/2)+1;
c_idx = round(eh.samples/2)+1;
imgAtSmallestWav = absImages(:,:,end);
[r,c,~] = ind2sub(size(imgAtSmallestWav),find(imgAtSmallestWav == max(imgAtSmallestWav(:))));

if isempty(x_pos) && isempty(y_pos)
    if add==1
        for i=1:numWav
            spec(i) = sum(sum(absImages(:,:,i)));
        end
    elseif add==2
        for i=1:numWav
            temp = absImages(:,:,i);
            spec(i) = temp(r_idx, c_idx);
        end
    else
        spec = squeeze(absImages(r,c,:));
        
    end
else
    spec = squeeze(absImages(x_pos,y_pos,:));
end
