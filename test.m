clear all;
%%
femmfunc(40)

   cfg = struct('nr',8,'nt',24,'r_inner',31.5,'r_outer',55,...
     'theta_span_deg',15,'slot_center_deg',7.5,...
     'lip_w',1.8,'lip_h',2.5,'top_w',8.7,'top_h',2,...
     'bottom_w',15.5,'bottom_h',12.645+2.5,'slot_depth',17.5,...
     'yoke_buffer_deg',0.3);
   domain = stator_design_domain(cfg);
   visualize_design_domain(domain); % 可视化确认掩膜是否符合预期
%%
% [T_avg, T_ripple] = scantorque(3.5)

% best_param = 25;  % = 39.3212
% 
% % 1. 看看最优点的指标
% femmfunc(best_param)
% [T_avg_best, T_ripple_best] = scantorque(3.5);
% F_best = T_avg_best - 1*T_ripple_best;