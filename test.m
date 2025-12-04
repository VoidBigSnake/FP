clear all;
%%
femmfunc(30)
%%
[T_avg, T_ripple] = scantorque(3.5)

% best_param = 25;  % = 39.3212
% 
% % 1. 看看最优点的指标
% femmfunc(best_param)
% [T_avg_best, T_ripple_best] = scantorque(3.5);
% F_best = T_avg_best - 1*T_ripple_best;