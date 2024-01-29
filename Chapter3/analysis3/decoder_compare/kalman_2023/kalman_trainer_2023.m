function [C,A,R,Q,P] = kalman_trainer_2023(xtrain,ztrainfilt,W)
% 26/05/2023

% adapted from kalman_trainer
% input ztrainfilt, assume filtered already

[n,train_length] = size(xtrain);

xtrain = xtrain';
C = eye(n);
X1 = xtrain(1:end-1,:);
X2 = xtrain(2:end,:);
A = X2'*X1*inv(X1'*X1);
R = (1/(train_length-1))*(X2'-A*X1')*(X2'-A*X1')';
Q = (1/train_length)*(ztrainfilt*W-xtrain*C)'*(ztrainfilt*W-xtrain*C);
P = eye(n);