function domain = stator_design_domain(cfg)
%STATORDESIGNDOMAIN 构建 15° 定子扇区的设计掩膜。
%  DOMAIN = STATORDESIGNDOMAIN(CFG) 返回网格边界与逻辑掩膜，
%  描述优化器可以编辑的极坐标单元。
%
%  CFG 必填字段：
%    nr             - 径向单元数（例：8）
%    nt             - 15° 内的切向单元数（例：20）
%    r_inner        - 定子内半径（气隙侧）
%    r_outer        - 定子外半径（轭外圆）
%    slot_r_inner   - 槽口起始半径（槽开口内缘）
%    slot_r_outer   - 槽口结束半径（槽开口外缘/线圈窗口）
%    slot_span_deg  - 15° 扇区内槽窗口的角度跨度
%    theta_span_deg - 设计扇区总角度（15）
%
%  CFG 可选字段：
%    coil_keepout_deg - 槽口两侧额外冻结的半角（默认 0）
%    yoke_buffer_deg  - 扇区边界附近的冻结半角，利于周期边界（默认 0）
%
%  输出 DOMAIN 字段：
%    r_edges, theta_edges - 网格边界
%    design_mask          - nt×nr 逻辑矩阵，true 表示可由优化器切换铁/空
%                            ；false 表示冻结单元
%    slot_mask            - nt×nr 逻辑矩阵，true 表示必须保持不变的槽/线圈区
%    notes                - 文字提醒的单元格数组
%
%  目标是在保持槽窗口与线圈区域不变的前提下，允许优化器修改铁轭与齿肩。
%
%  示例：
%    cfg = struct('nr',8,'nt',20,'r_inner',22,'r_outer',40,...
%                 'slot_r_inner',22.5,'slot_r_outer',28,...
%                 'slot_span_deg',6,'theta_span_deg',15,...
%                 'coil_keepout_deg',1.0,'yoke_buffer_deg',0.5);
%    domain = stator_design_domain(cfg);
%    visualize_design_domain(domain); % 可视化确认掩膜是否符合预期
%
%  参见 APPLY_DESIGN_MASK。

% 使用较老版本 MATLAB/Octave 时，避免 ARGUMENTS 语法带来的解析错误，改为
% 传统的字段检查与默认值填充。
required_fields = { ...
    'nr', 'nt', 'r_inner', 'r_outer', ...
    'slot_r_inner', 'slot_r_outer', 'slot_span_deg', 'theta_span_deg'};
for f = required_fields
    if ~isfield(cfg, f{1})
        error('缺少必填字段 cfg.%s', f{1});
    end
end

% 默认值（如已提供则不覆盖）。
if ~isfield(cfg, 'coil_keepout_deg'); cfg.coil_keepout_deg = 0; end
if ~isfield(cfg, 'yoke_buffer_deg');  cfg.yoke_buffer_deg  = 0; end

% 数值检查（正数/非负）。
mustBePositive = @(v,name) assert(isnumeric(v) && isscalar(v) && v>0, ...
    '字段 %s 需为正数标量，当前为 %s', name, mat2str(v));
mustBeNonnegative = @(v,name) assert(isnumeric(v) && isscalar(v) && v>=0, ...
    '字段 %s 需为非负标量，当前为 %s', name, mat2str(v));

mustBePositive(cfg.nr, 'nr');
mustBePositive(cfg.nt, 'nt');
mustBePositive(cfg.r_inner, 'r_inner');
mustBePositive(cfg.r_outer, 'r_outer');
mustBePositive(cfg.slot_r_inner, 'slot_r_inner');
mustBePositive(cfg.slot_r_outer, 'slot_r_outer');
mustBePositive(cfg.slot_span_deg, 'slot_span_deg');
mustBePositive(cfg.theta_span_deg, 'theta_span_deg');
mustBeNonnegative(cfg.coil_keepout_deg, 'coil_keepout_deg');
mustBeNonnegative(cfg.yoke_buffer_deg, 'yoke_buffer_deg');

if cfg.r_outer <= cfg.r_inner
    error('r_outer must be greater than r_inner (got %.2f vs %.2f).', ...
          cfg.r_outer, cfg.r_inner);
end

if cfg.slot_r_outer <= cfg.slot_r_inner
    error('slot_r_outer must be greater than slot_r_inner (got %.2f vs %.2f).', ...
          cfg.slot_r_outer, cfg.slot_r_inner);
end

% 构建网格边界与中心。
r_edges = linspace(cfg.r_inner, cfg.r_outer, cfg.nr + 1);
theta_edges = linspace(0, cfg.theta_span_deg, cfg.nt + 1); % 角度制
[theta_centers, r_centers] = meshgrid(
    theta_edges(1:end-1) + diff(theta_edges)/2,
    r_edges(1:end-1) + diff(r_edges)/2);
% theta_centers 为 nt×1，并在行方向复制；转置后方便掩膜运算。
theta_centers = theta_centers';
r_centers = r_centers';

% 槽窗口掩膜（围绕扇区中心的切向跨度）。
slot_half = cfg.slot_span_deg/2 + cfg.coil_keepout_deg;
slot_center = cfg.theta_span_deg/2;
theta_slot = abs(theta_centers - slot_center) <= slot_half;
r_slot = (r_centers >= cfg.slot_r_inner) & (r_centers <= cfg.slot_r_outer);
slot_mask = theta_slot & r_slot;

% 扇区边界附近的缓冲区，保持周期边界节点稳定。
edge_mask = (theta_centers <= cfg.yoke_buffer_deg) | ...
            (theta_centers >= (cfg.theta_span_deg - cfg.yoke_buffer_deg));

% 设计掩膜：仅允许在轭/齿肩编辑，不包含槽与边界缓冲区。
design_mask = ~slot_mask & ~edge_mask;

% 汇总结果。
domain = struct();
domain.r_edges = r_edges;
domain.theta_edges = theta_edges;
domain.design_mask = design_mask;
domain.slot_mask = slot_mask;
domain.notes = {
    'design_mask: true = 优化器可切换铁/空气，false = 冻结';
    'slot_mask: true = 槽/线圈区保持不变';
    'apply_design_mask(bits, domain, base_val) 会强制冻结区域';
    sprintf('sector: %.1f deg, slot span: %.1f deg (keepout %.1f deg)', ...
            cfg.theta_span_deg, cfg.slot_span_deg, cfg.coil_keepout_deg)
};

end

function grid = apply_design_mask(bits, domain, base_val)
%APPLY_DESIGN_MASK 将基因比特映射到网格并冻结槽/缓冲区。
%  GRID = APPLY_DESIGN_MASK(BITS, DOMAIN, BASE_VAL) 返回 nt×nr 数组，
%  可编辑单元取自 BITS（重塑后），冻结区域（槽/线圈或边界缓冲）
%  则被强制为 BASE_VAL（例如铁 = 1）。
%
%  BITS 需为长度 nt*nr 的向量，与 DOMAIN 掩膜尺寸一致。

if numel(bits) ~= numel(domain.design_mask)
    error('Bitstring length %d does not match grid size %d.', ...
          numel(bits), numel(domain.design_mask));
end

grid = reshape(bits, size(domain.design_mask));
freeze_mask = ~domain.design_mask; % 包含槽与边界缓冲
grid(freeze_mask) = base_val;
end
