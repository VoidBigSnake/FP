function visualize_design_domain(domain)
%VISUALIZE_DESIGN_DOMAIN 可视化 15° 扇区的设计掩膜与冻结区域。
%  VISUALIZE_DESIGN_DOMAIN(DOMAIN) 使用 imagesc 显示 design_mask（可编辑）
%  与 slot_mask（槽/线圈冻结区），帮助确认 stator_design_domain(cfg)
%  生成的掩膜是否生效。
%
%  颜色说明：
%    绿色 = 设计域（true in design_mask）
%    灰色 = 槽/线圈或边界缓冲（slot_mask 或 edge 区域）
%
%  提示：掩膜仅用于优化阶段的网格/基因映射，并不会直接修改 FEMM
%  几何；若希望在 FEMM 中看到结构变化，需要在生成几何时按掩膜
%  写入或移除铁/空气区域。

if ~isstruct(domain) || ~isfield(domain,'design_mask') || ~isfield(domain,'slot_mask')
    error('输入需为包含 design_mask 和 slot_mask 的结构体。');
end

mask = domain.design_mask;
slot = domain.slot_mask;

% 构造双通道掩膜：优先显示槽/冻结区，其余为设计域。
plot_mask = zeros(size(mask));
plot_mask(mask) = 1;      % 设计域
plot_mask(slot) = -1;     % 槽/线圈冻结区

% 以物理尺度显示，避免索引坐标导致的拉伸：
%   x 轴 = 切向弧长（以内半径为基准），y 轴 = 径向距离。
r_edges = domain.r_edges;
theta_edges = domain.theta_edges;
r_base = r_edges(1);

u_min = deg2rad(theta_edges(1)) * r_base;
u_max = deg2rad(theta_edges(end)) * r_base;
v_min = 0;
v_max = r_edges(end) - r_base;

imagesc([u_min, u_max], [v_min, v_max], plot_mask');
axis equal tight; axis ij;
colormap([0.5 0.5 0.5; 0 0.7 0]);
colorbar('Ticks',[-1 1],'TickLabels',{'冻结区','设计域'});
title('设计域掩膜（15° 扇区）');
xlabel('切向弧长 / mm（r = r_{in}）');
ylabel('径向距离 / mm');

% 如提供槽轮廓（局部 u/v），在图上覆盖轮廓线，便于对齐检查。
if isfield(domain, 'slot_outline_uv')
    hold on;
    u = domain.slot_outline_uv.u;
    v = domain.slot_outline_uv.v;
    if isfield(domain, 'slot_center_deg')
        slot_center_deg = domain.slot_center_deg;
    else
        slot_center_deg = (theta_edges(1) + theta_edges(end)) / 2;
    end

    % 局部坐标 -> 物理弧长/径向距离，与 imagesc 的坐标一致。
    r = r_base + v;
    theta = slot_center_deg + rad2deg(u ./ max(r, eps));

    u_abs = deg2rad(theta) * r_base; % 以内半径为基准的弧长
    v_abs = v;

    plot(u_abs, v_abs, 'k-', 'LineWidth', 1.2, 'DisplayName', '槽轮廓参考');
    hold off;
end

gridsize = size(mask);
text(u_min, v_max + (v_max - v_min)*0.02, sprintf('nt=%d, nr=%d', gridsize(1), gridsize(2)));

end