function [Et_ff, Ef_ff] = getFarField(Et_bpf, Ef_bpf,startIdx, endIdx)

%get far field
%first crop the filtered near-field image of the source and scattered fields
Et_ff=Et_bpf(startIdx:endIdx,startIdx:endIdx);
Et_ff(isnan(Et_ff))=0;     Et_ff(isinf(Et_ff))=0;
Ef_ff= Ef_bpf(startIdx:endIdx, startIdx:endIdx);
Ef_ff(isnan(Ef_ff))=0;     Ef_ff(isinf(Ef_ff))=0;
