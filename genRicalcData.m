function genRicalcData(hdrPath, absImgPath, outputPath)
% generate file for ricalc
[wav, absImages, spec] = loadFTIR(hdrPath, absImgPath, 3);
wavAbs = [wav' spec];
csvwrite(outputPath,wavAbs);