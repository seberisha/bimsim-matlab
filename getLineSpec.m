function lsp = getLineSpec(i)
lspec = {'-','--',':','-.'};
lsp = lspec{mod(i,numel(lspec))+1};