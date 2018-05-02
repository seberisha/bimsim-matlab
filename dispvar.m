function dispvar(varName, varValue)
% fprintf('%s = ',inputname(1)) ;
% disp(var) ;
if nargin==0
    fprintf('\n ******* \n\n')
elseif nargin==2
    fprintf('%s = %g\n', varName, varValue);
end
%  disp(sprintf('%s = %g\n', varName, varValue));
