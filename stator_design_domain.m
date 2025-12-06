function domain = stator_design_domain(cfg)
%STATORDESIGNDOMAIN 构建 15° 定子扇区的设计掩膜。
%  DOMAIN = STATORDESIGNDOMAIN(CFG) 返回网格边界与逻辑掩膜，
%  描述优化器可以编辑的极坐标单元。
%
%  两种模式：
%  1) 'slot'（默认）：槽口矩形 + 六边形槽体被冻结，只在铁轭/齿肩内优化。
%     适用于“槽不动”的拓扑优化。
%  2) 'yoke'：直接用径向内外半径限定 15° 扇区的铁轭环带，忽略槽型细节，
%     适合“只在铁轭”里开窗口或挖孔的场景；如希望自动跳过槽/线圈，可额外
%     提供槽几何字段（与 slot 模式相同）以自动冻结槽腔，这样 design_mask
%     仅覆盖铁轭区域。
%
%  公共必填字段：
%    nr, nt                - 径向/切向单元数（例：8、20）
%    r_inner, r_outer      - 定子内/外半径（mm）
%    theta_span_deg        - 设计扇区总角度（度，例：15）
%
%  模式 'slot' 额外必填字段：
%    slot_center_deg       - 槽中心角（度，通常=theta_span_deg/2）
%    mouth_w, mouth_h      - 槽口矩形的宽/高（mm）
%    top_w,  top_h         - 槽体上边宽、上边到槽口顶的高（mm）
%    bottom_w, bottom_h    - 槽体底边宽、槽体高度（mm）
%    slot_depth            - 槽总深度（mm，含底部外延）
%    lip_w                 - 齿唇宽（mm，对应槽口两侧微收口）
%    bottom_ext            - 槽底向外额外外延量（mm，结构图中 0.5）
%
%  模式 'yoke' 额外必填字段：
%    design_r_inner        - 设计域内径（mm），应为铁轭开始半径
%    design_r_outer        - 设计域外径（mm），不超过 r_outer
%    （可选）若同时给出槽几何字段与 slot_center_deg，将自动冻结槽腔
%
%  CFG 可选字段：
%    yoke_buffer_deg       - 扇区边界附近的冻结半角，利于周期边界（默认 0）
%    mode                  - 'slot' 或 'yoke'（默认 'slot'）
%
%  输出 DOMAIN 字段：
%    r_edges, theta_edges  - 网格边界
%    design_mask           - nt×nr 逻辑矩阵，true = 可由优化器切换铁/空；false = 冻结
%    slot_mask             - nt×nr 逻辑矩阵，true = 槽/线圈/冻结区
%    slot_outline_uv       - 槽截面顶点（仅 slot 模式下给出，便于调试/复用）
%    notes                 - 文字提醒
%
%  示例（槽冻结，优化轭/齿肩）：
%    cfg = struct('nr',8,'nt',24,'r_inner',31.5,'r_outer',55,...
%      'theta_span_deg',15,'slot_center_deg',7.5,...
%      'mouth_w',8.7,'mouth_h',2,'top_w',8.7,'top_h',2,...
%      'bottom_w',15.5,'bottom_h',10.145,'slot_depth',17.5,...
%      'lip_w',1.8,'bottom_ext',0.5,'yoke_buffer_deg',0.3);
%    domain = stator_design_domain(cfg);
%    visualize_design_domain(domain); % 可视化确认掩膜是否符合预期
%
%  示例（只在铁轭环带内优化，忽略槽形）：
%    cfg = struct('nr',8,'nt',24,'r_inner',31.5,'r_outer',55,...
%      'theta_span_deg',15,'design_r_inner',48,'design_r_outer',55,...
%      'mode','yoke','yoke_buffer_deg',0.3);
%    domain = stator_design_domain(cfg);
%    visualize_design_domain(domain);
%
%  示例（在铁轭环带内优化，同时自动识别槽并冻结）：
%    cfg = struct('nr',8,'nt',24,'r_inner',31.5,'r_outer',55,...
%      'theta_span_deg',15,'design_r_inner',48,'design_r_outer',55,...
%      'mode','yoke','yoke_buffer_deg',0.3,'slot_center_deg',7.5,...
%      'mouth_w',8.7,'mouth_h',2,'top_w',8.7,'top_h',2,...
%      'bottom_w',15.5,'bottom_h',10.145,'slot_depth',17.5,...
%      'lip_w',1.8,'bottom_ext',0.5);
%    domain = stator_design_domain(cfg);
%    visualize_design_domain(domain);
%
%  参见 APPLY_DESIGN_MASK、VISUALIZE_DESIGN_DOMAIN。

% 使用较老版本 MATLAB/Octave 时，避免 ARGUMENTS 语法带来的解析错误，改为
% 传统的字段检查与默认值填充。
if ~isfield(cfg, 'mode'); cfg.mode = 'slot'; end
if ~isfield(cfg, 'yoke_buffer_deg');  cfg.yoke_buffer_deg  = 0; end

% 必填字段按模式拆分。
common_required = {'nr','nt','r_inner','r_outer','theta_span_deg'};
slot_required = { ...
    'slot_center_deg','mouth_w','mouth_h','top_w','top_h', ...
    'bottom_w','bottom_h','slot_depth','lip_w','bottom_ext'};
yoke_required = {'design_r_inner','design_r_outer'};

for f = common_required
    if ~isfield(cfg, f{1})
        error('缺少必填字段 cfg.%s', f{1});
    end
end

switch lower(cfg.mode)
    case 'slot'
        for f = slot_required
            if ~isfield(cfg, f{1})
                error('模式 slot 需提供 cfg.%s', f{1});
            end
        end
    case 'yoke'
        for f = yoke_required
            if ~isfield(cfg, f{1})
                error('模式 yoke 需提供 cfg.%s', f{1});
            end
        end
    otherwise
        error('未知模式 cfg.mode=%s（支持 ''slot'' 或 ''yoke''）', cfg.mode);
end

% 默认值不覆盖用户输入。

% 数值检查（正数/非负）。
mustBePositive = @(v,name) assert(isnumeric(v) && isscalar(v) && v>0, ...
    '字段 %s 需为正数标量，当前为 %s', name, mat2str(v));
mustBeNonnegative = @(v,name) assert(isnumeric(v) && isscalar(v) && v>=0, ...
    '字段 %s 需为非负标量，当前为 %s', name, mat2str(v));

mustBePositive(cfg.nr, 'nr');
mustBePositive(cfg.nt, 'nt');
mustBePositive(cfg.r_inner, 'r_inner');
mustBePositive(cfg.r_outer, 'r_outer');
mustBePositive(cfg.theta_span_deg, 'theta_span_deg');
mustBeNonnegative(cfg.yoke_buffer_deg, 'yoke_buffer_deg');

if cfg.r_outer <= cfg.r_inner
    error('r_outer must be greater than r_inner (got %.2f vs %.2f).', ...
          cfg.r_outer, cfg.r_inner);
end

if strcmpi(cfg.mode, 'slot')
    mustBePositive(cfg.mouth_w, 'mouth_w');
    mustBePositive(cfg.mouth_h, 'mouth_h');
    mustBePositive(cfg.top_w, 'top_w');
    mustBePositive(cfg.top_h, 'top_h');
    mustBePositive(cfg.bottom_w, 'bottom_w');
    mustBePositive(cfg.bottom_h, 'bottom_h');
    mustBePositive(cfg.slot_depth, 'slot_depth');
    mustBePositive(cfg.lip_w, 'lip_w');
    mustBePositive(cfg.bottom_ext, 'bottom_ext');

    if cfg.slot_depth <= cfg.mouth_h
        error('slot_depth must exceed mouth_h (got %.2f vs %.2f).', ...
              cfg.slot_depth, cfg.mouth_h);
    end

    if cfg.slot_depth <= cfg.top_h + cfg.bottom_h
        error('slot_depth 应大于 top_h+bottom_h，当前 %.2f <= %.2f', ...
              cfg.slot_depth, cfg.top_h + cfg.bottom_h);
    end
else
    mustBePositive(cfg.design_r_inner, 'design_r_inner');
    mustBePositive(cfg.design_r_outer, 'design_r_outer');

    if cfg.design_r_outer > cfg.r_outer
        error('design_r_outer=%.2f 不应超过 r_outer=%.2f', ...
              cfg.design_r_outer, cfg.r_outer);
    end

    if cfg.design_r_inner < cfg.r_inner
        error('design_r_inner=%.2f 不应小于 r_inner=%.2f', ...
              cfg.design_r_inner, cfg.r_inner);
    end

    if cfg.design_r_outer <= cfg.design_r_inner
        error('design_r_outer 必须大于 design_r_inner (%.2f vs %.2f)', ...
              cfg.design_r_outer, cfg.design_r_inner);
    end
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

% 扇区边界附近的缓冲区，保持周期边界节点稳定。
edge_mask = (theta_centers <= cfg.yoke_buffer_deg) | ...
            (theta_centers >= (cfg.theta_span_deg - cfg.yoke_buffer_deg));

slot_outline = struct('u',[],'v',[]);
slot_mask = false(size(theta_centers));
slot_frozen_in_yoke = false;

if strcmpi(cfg.mode, 'slot')
    % 槽截面（局部 u/v，多边形顶点顺时针）。
    slot_outline = build_slot_outline(cfg);

    % 将极坐标网格转换为槽局部 u/v：
    % u = r*dtheta（切向弧长），v = r - r_inner（径向到齿根距离）。
    dtheta_rad = deg2rad(theta_centers - cfg.slot_center_deg);
    u_local = r_centers .* dtheta_rad;
    v_local = r_centers - cfg.r_inner;

    % 槽/线圈冻结区：在槽口矩形 + 六边形槽体多边形内部。
    in_slot = inpolygon(u_local(:), v_local(:), slot_outline.u, slot_outline.v);
    slot_mask = reshape(in_slot, size(theta_centers));

    % 设计掩膜：仅允许在轭/齿肩编辑，不包含槽与边界缓冲区。
    design_mask = ~slot_mask & ~edge_mask;
else
    % Yoke 模式：仅在指定的径向环带内可编辑，其余保持冻结。
    radial_band = (r_centers >= cfg.design_r_inner) & ...
                  (r_centers <= cfg.design_r_outer);
    design_mask = radial_band & ~edge_mask;

    % 可选：如果提供了完整槽几何字段，则在 yoke 模式下也自动冻结槽腔，
    % 确保 design_mask 只作用在铁轭而不触碰槽/线圈。
    slot_shape_fields = {
        'mouth_w','mouth_h','top_w','top_h',...
        'bottom_w','bottom_h','slot_depth','lip_w','bottom_ext'
    };
    has_slot_shape = all(cellfun(@(f)isfield(cfg,f), slot_shape_fields));
    if has_slot_shape
        if ~isfield(cfg, 'slot_center_deg')
            slot_center = cfg.theta_span_deg/2; % 缺省取扇区中心
        else
            slot_center = cfg.slot_center_deg;
        end

        tmp_cfg = cfg; % 不修改原结构
        tmp_cfg.slot_center_deg = slot_center;
        slot_outline = build_slot_outline(tmp_cfg);

        dtheta_rad = deg2rad(theta_centers - slot_center);
        u_local = r_centers .* dtheta_rad;
        v_local = r_centers - cfg.r_inner;
        in_slot = inpolygon(u_local(:), v_local(:), slot_outline.u, slot_outline.v);
        slot_mask = reshape(in_slot, size(theta_centers));
        design_mask = design_mask & ~slot_mask;
        slot_frozen_in_yoke = true;
    else
        % 将环带外部视为“冻结区”，便于可视化。
        slot_mask = ~radial_band | edge_mask;
    end
end

% 汇总结果。
domain = struct();
domain.r_edges = r_edges;
domain.theta_edges = theta_edges;
domain.design_mask = design_mask;
domain.slot_mask = slot_mask;
domain.slot_outline_uv = slot_outline;
domain.notes = {
    'design_mask: true = 优化器可切换铁/空气，false = 冻结';
    'slot_mask: true = 槽/线圈或非设计区保持不变';
    sprintf('模式=%s, 扇区 %.1f°', cfg.mode, cfg.theta_span_deg)
};
if strcmpi(cfg.mode, 'yoke')
    if slot_frozen_in_yoke
        domain.notes{end+1} = 'yoke 模式已按槽几何自动冻结槽腔，仅铁轭可编辑';
    else
        domain.notes{end+1} = 'yoke 模式未提供槽几何，未额外识别槽腔';
    end
end

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

function outline = build_slot_outline(cfg)
%BUILD_SLOT_OUTLINE 生成槽口矩形 + 六边形槽体的局部坐标（u/v）。
%  基于 strukturFemm 的槽型：
%    - mouth_w / mouth_h     : 槽口矩形
%    - top_w  / top_h        : 槽体上边宽与槽口顶到上边的高
%    - bottom_w / bottom_h   : 槽体底边宽与槽体高度
%    - slot_depth            : 槽总深度（含 bottom_ext）
%    - bottom_ext (0.5 mm)   : 底部外延
%    - lip_w                 : 齿唇宽，控制槽口两侧向内收一点
%
%  坐标系：v 从齿根圆向外（mm），u 为切向位移（mm，右正）。

v0 = 0;                        % 齿根处
v_mouth = v0 + cfg.mouth_h;    % 槽口矩形上边
v_top   = v_mouth + cfg.top_h; % 槽体上边高度
v_bottom= v_top + cfg.bottom_h;% 槽体底边高度
v_tip   = cfg.slot_depth + cfg.bottom_ext; % 槽底最外侧

w_mouth = cfg.mouth_w/2;
w_top   = cfg.top_w/2;
w_bottom= cfg.bottom_w/2;
w_lip   = cfg.lip_w/2; % 齿唇微收口（v=0 处）

outline_u = [ ...
    w_lip, v0;          % 齿唇处（略窄）
    w_mouth, v_mouth;   % 槽口矩形右上
    w_top, v_top;       % 槽体上边右端
    w_bottom, v_bottom; % 槽体底边右端
    0, v_tip;           % 底尖
   -w_bottom, v_bottom; % 底边左端
   -w_top, v_top;       % 上边左端
   -w_mouth, v_mouth;   % 槽口矩形左上
   -w_lip, v0           % 齿唇处（略窄）
];

% 变为行向量，便于 inpolygon 调用。
outline.u = outline_u(:,1)';
outline.v = outline_u(:,2)';
end
