clear all;
% %%
% femmfunc_design();
% 
%    cfg = struct('nr',8,'nt',16,'r_inner',31.5,'r_outer',54,...
%      'theta_span_deg',15,'slot_center_deg',7.5,...
%      'lip_w',1.8,'lip_h',2.5,'top_w',8.7,'top_h',2,...
%      'bottom_w',15.5,'bottom_h',12.645+2.5,'slot_depth',17.5,...
%      'yoke_buffer_deg',0.1);
%    domain0 = stator_design_domain_blank(cfg);
%    % visualize_design_domain(domain); % 可视化确认掩膜是否符合预期
%    theta_offset_deg = 0; 
%    % femm_debug_design_domain(domain, cfg, theta_offset_deg);
%    femm_draw_design_grid(domain0, theta_offset_deg, 20);
% % visualize_design_cells(domain, theta_offset_deg);
% %%
% domain = prepare_design_domain(domain0, theta_offset_deg);
% Nd      = domain.Nd;         % 可设计格子数
% r_c     = hypot(domain.x_c, domain.y_c);   % 每个格子中心的半径
% circNames = {'A+','A-','B+','B-','C+','C-'};  % 6 条电路
% % sectorCircuit(s) = 第 s 个 15° 扇区的电路编号（对应 circNames）
% % 例如：0~15° → A+, 15~30° → B+, 30~45° → C+, 45~60° → A-, 60~75° → B-, 75~90° → C-
% sectorCircuit = [4 3 6 5 2 1];    % 按你自己的相序填；这里只是一个例子
% 
% phase_id_sector = cell(6,1);      % 6 个扇区，每个一个 Nd×1 向量
% 
% for s = 1:6
%     pid = sectorCircuit(s);       % 1..6
%     phase_id_sector{s} = pid * ones(Nd,1);   % 整个扇区里的铜都挂到同一条电路
% end
% 
% Nd   = domain.Nd;
% L    = 2 * Nd;
% 
% % 1) 随机一条基因串
% bits = rand(1, L) > 0.5;
% 
% % 2) 6 个扇区的旋转角
% sector_offsets_deg = 0:15:75;
% 
% % 3) 材料名（按你 FEMM 工程改）
% mats = struct();
% mats.air      = 'Air';
% mats.iron     = 'Pure Iron';
% mats.copper   = 'Copper';
% mats.meshSize = 1.0;
% 
% femm_apply_design_bits_rep6(bits, domain, sector_offsets_deg, ...
%                             phase_id_sector, mats, circNames, 10);

%%
% [T_avg, T_ripple] = scantorque(3.5)

% best_param = 25;  % = 39.3212
% 
% % 1. 看看最优点的指标
% femmfunc(best_param)
% [T_avg_best, T_ripple_best] = scantorque(3.5);
% F_best = T_avg_best - 1*T_ripple_best;

%%
% b1 = randi([0,1], 1, 2*domain.Nd);
% b2 = randi([0,1], 1, 2*domain.Nd);
% 
% [J1,Tavg1,Trip1] = eval_design_femm(b1, ctx);
% [J2,Tavg2,Trip2] = eval_design_femm(b2, ctx);
% 
% disp([J1 Tavg1 Trip1]);
% disp([J2 Tavg2 Trip2]);
%% ---------- 0. 一次性准备设计域 + ctx ----------
clear; clc;

% ---- 0.1 设计域参数（示意，改成你自己的 cfg）----
   cfg = struct('nr',9,'nt',5,'r_inner',31.5,'r_outer',52,...
     'theta_span_deg',15,'slot_center_deg',7.5,...
     'lip_w',1.8,'lip_h',2.5,'top_w',8.7,'top_h',2,...
     'bottom_w',15.5,'bottom_h',12.645+2.5,'slot_depth',17.5,...
     'yoke_buffer_deg',0.1);

domain0 = stator_design_domain_blank(cfg);  % 你现有的"空白定子设计域"函数
theta0_deg = 0;                             % 第一块 15° 就是 0~15°
domain  = prepare_design_domain(domain0, theta0_deg);  % 加上 x_c,y_c,Nd 等

Nd   = domain.Nd;
Lbit = Nd;

% ---- 0.2 6 个扇区的相别模板（你应该已经有 build_phase_id_sector）----
circNames = {'A+','A-','B+','B-','C+','C-'};  % 6 条电路
sectorCircuit = [4 3 6 5 2 1];    % 按你自己的相序填；这里只是一个例子
phase_id_sector = cell(6,1);      % 6 个扇区，每个一个 Nd×1 向量

for s = 1:6
    pid = sectorCircuit(s);       % 1..6
    phase_id_sector{s} = pid * ones(Nd,1);   % 整个扇区里的铜都挂到同一条电路
end

% ---- 0.3 材料 & 电路名 ----
mats = struct();
mats.air      = 'Air';
mats.iron     = 'Pure Iron';
mats.copper   = 'Copper';
mats.meshSize = 1.0;   % 设计域网格尺寸，自己调

% ---- 0.4 每个电路的"总匝数"设定（例子：每个电路 100 匝）----
N_phase_total = 100 * ones(1, numel(circNames));

% ---- 0.5 FEMM 模板文件 & 分组号 ----
ctx = struct();

ctx.domain             = domain;
ctx.sector_offsets_deg = 0:15:75;        % 或者 30+0:15:75 之类
ctx.phase_id_sector    = phase_id_sector;
ctx.groupId_design     = 30;             % 设计域 group 号你自己定

ctx.baseFemFile = 'blank_9x5.fem';

ctx.mats = struct();
ctx.mats.air      = 'Air';
ctx.mats.iron     = 'Pure Iron';
ctx.mats.copper   = 'Copper';
ctx.mats.meshSize = 4.0;

ctx.circNames = {'A+','A-','B+','B-','C+','C-'};

ctx.T_min       = 0.5;    % 目标平均转矩要求（先随便放一个）
ctx.penaltyCoef = 10;

ctx.N_phase_total   = N_phase_total;

ctx.groupId_core = 30;
ctx.groupId_ring = 31;
ctx.inset_r_ratio  = 0.2;
ctx.inset_th_ratio = 0.2;

ctx.cfg = cfg;

% 你 trit 的材料编码（按你实际改）
ctx.airCode  = 0;
ctx.ironCode = 1;
ctx.cuCode   = 2;

% 你的拓扑是扇区重复的，所以建议 θ 方向周期连通
ctx.thetaPeriodic = true;


%% ---------- 1. 用一条 bits 测试 eval_design_femm ----------
Nd   = domain.Nd;
Lbit = Nd;

S = load('seed_bits_0_15deg_1.mat', 'seed_bits', 'cfg');
assert(numel(S.seed_bits) == ctx.domain.Nd, 'seed_bits length mismatch: cfg changed?');
seed_bits = S.seed_bits;

[J, T_avg, T_ripple] = eval_design_femm(seed_bits, ctx);

fprintf('Test design: J = %.4f, T_avg = %.3f Nm, T_ripple = %.3f Nm\n', ...
        J, T_avg, T_ripple);

%%
refFemFile = 'best_J0.150072_20251230_212929.fem';
openfemm;
opendocument(refFemFile)
[T_avg, T_ripple] = scantorque(3.5);

