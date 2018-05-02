%%Gaussian Maximum Likelihood Classifier
%%Case 3, arbitrary covariance matrix
 
function [Y_Test] = GaussianMLClassifier3(X_Train, Y_Train, X_Test)
    ClassNum = 3;                         
    d = size(X_Train,2);  
    for k = 1 : ClassNum                   %Train Mu and Sigma for each class
        index = find( Y_Train == k);       %By using funtion of HW1
        [Mu(k,:), Sigma(:,:,k)] = GaussianMLEstimator(X_Train(index,:));
    end
    for j = 1 : length(X_Test)             %for each sample of X_Test,
        for i = 1 : ClassNum               %get its value of discriminant funtion of each class
            G_part1 = (X_Test(j,:) - Mu(i,:))*Sigma(:,:,i)^(-1)*(X_Test(j,:) - Mu(i,:))';
            G_part2 = -0.5*d*log(2*pi) - 0.5*log(det(Sigma(:,:,i))) + log(1/ClassNum);
            G(i) = -0.5*G_part1 + G_part2;
        end
        [~, b] = max(G);                   %get the maximum value of g(i)
        Y_Test(j) = b;                     %i is the class number this sample belongs to
    end
