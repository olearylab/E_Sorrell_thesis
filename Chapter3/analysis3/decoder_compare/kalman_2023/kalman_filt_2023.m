function [xpost,P] = kalman_filt_2023(C,A,R,Q,P,xpred,xpost,gain_weight)

xpost = A*xpost;
P = A*P*A' + R;
K = gain_weight*P*C*inv(C'*P*C + Q);
xpost = xpost + K*(xpred' - C*xpost);
P = (eye(length(xpost)) - K*C)*P;

