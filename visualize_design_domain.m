function visualize_design_domain(domain)
%VISUALIZE_DESIGN_DOMAIN 可视化 15° 扇区的设计掩膜与冻结区域。
%  VISUALIZE_DESIGN_DOMAIN(DOMAIN) 使用 imagesc 显示 design_mask（可编辑）
%  与 slot_mask（槽/线圈或非设计区），帮助确认 stator_design_domain(cfg)
%  生成的掩膜是否生效（支持 slot / yoke 两种模式）。
%
%  颜色说明：
%    绿色 = 设计域（true in design_mask）
%    灰色 = 槽/线圈或边界缓冲/非设计环带（slot_mask 或 edge 区域）
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

imagesc(plot_mask');
axis equal tight; axis ij;
colormap([0.5 0.5 0.5; 0 0.7 0]);
colorbar('Ticks',[-1 1],'TickLabels',{'冻结区','设计域'});
title('设计域掩膜（15° 扇区）');
xlabel('切向单元索引');
ylabel('径向单元索引');

% 如提供槽轮廓（局部 u/v），在图上覆盖轮廓线，便于对齐检查。
if isfield(domain, 'slot_outline_uv')
    hold on;
    u = domain.slot_outline_uv.u;
    v = domain.slot_outline_uv.v;
    plot((u - min(u)) / (max(u)-min(u)+eps) * (size(mask,1)-1) + 1, ...
         (v - min(v)) / (max(v)-min(v)+eps) * (size(mask,2)-1) + 1, ...
         'k-', 'LineWidth', 1.2, 'DisplayName', '槽轮廓参考');
    hold off;
end

gridsize = size(mask);
text(1, gridsize(2)+0.3, sprintf('nt=%d, nr=%d', gridsize(1), gridsize(2)));

end
