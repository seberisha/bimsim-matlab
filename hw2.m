load /media/sberisha/E4D3-804B/data.mat

ClassNum = 3;
FHalf1 = data.f1(1:500,:);
FHalf2 = data.f2(1:500,:);
FHalf3 = data.f3(1:500,:);                  %get the former half of dataset as
X_Train = cat(1,FHalf1,FHalf2,FHalf3);      %X training data
 
YHalf1 = zeros(500,1)+1;
YHalf2 = zeros(500,1)+2;
YHalf3 = zeros(500,1)+3;                    %making up the label value as
Y_Train = cat(1,YHalf1,YHalf2,YHalf3);      %Y training data
 
XHalf1 = data.f1(501:1000,:);
XHalf2 = data.f2(501:1000,:);
XHalf3 = data.f3(501:1000,:);               %get the former half of dataset as
X_Test = cat(1,XHalf1,XHalf2,XHalf3);       %X testing data
 
d = size(X_Train,2);
 
[Y_Test] = GaussianMLClassifier3(X_Train, Y_Train, X_Test);
 
Error = find( (Y_Test)' ~= Y_Train);
ErrorNum = length(Error);
Accuracy = 1 - ErrorNum/length(Y_Train);    %calculating accuracy