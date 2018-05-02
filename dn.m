function dn(varname)
  disp(sprintf('%s = %g\n', varname, evalin('caller',varname)));
end