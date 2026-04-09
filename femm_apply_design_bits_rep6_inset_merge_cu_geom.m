function femm_apply_design_bits_rep6_inset_merge_cu_geom( ...
    gene, domain, ctx, phase_id_sector, mats, circNames, ...
    groupIdCore, groupIdRing, groupIdCuGeom, turns_per_circ)

% 目标：先按网格连通关系得到铜岛，再仅对"铜岛外边界"做小幅内缩（默认 0.2 个格宽）
% - 铜岛内部共享边不内缩 => 不会出现你提到的内部缝隙
% - 铜岛外边界内缩 => 形成环绕铜岛的空气带

Nd = domain.Nd;
mat_code = decode_material_bits(gene, Nd);  % 0 air, 1 iron, 2 copper

alpha = 0.2;
if isfield(ctx, 'cu_island_inset_ratio')
    alpha = ctx.cu_island_inset_ratio;
elseif isfield(ctx, 'cu_inset_ratio')
    alpha = ctx.cu_inset_ratio;
elseif isfield(ctx, 'inset_r_ratio')
    alpha = ctx.inset_r_ratio;
end
alpha = max(0, min(alpha, 0.45));

% ===== 1) 删旧 label =====
mi_selectgroup(groupIdCore); mi_deleteselected();
mi_selectgroup(groupIdRing); mi_deleteselected();
mi_clearselected();

% ===== 2) 删旧动态铜几何 =====
mi_selectgroup(groupIdCuGeom);
mi_deleteselected();
mi_clearselected();

% ===== 基域点 =====
x0c = domain.x_core(:); y0c = domain.y_core(:);
x0r = domain.x_ring(:); y0r = domain.y_ring(:);
r0c = hypot(x0c,y0c); th0c = atan2d(y0c,x0c);
r0r = hypot(x0r,y0r);

c0 = mean(th0c);
d  = th0c - c0;

theta_edges = domain.theta_edges(:);
th_span = theta_edges(end) - theta_edges(1);
nSector = numel(phase_id_sector);

nr = ctx.cfg.nr; nt = ctx.cfg.nt;
assert(Nd == nr*nt, 'Nd != nr*nt');
r_edges = domain.r_edges(:);

% ===== k -> (ir,it) =====
th0c_wrap = wrap_to_span(th0c, theta_edges(1), theta_edges(end));
epsr = 1e-9; epsth = 1e-9;
r0c2 = r0c; r0c2(r0c2 >= r_edges(end)) = r_edges(end) - epsr;
th0c2 = th0c_wrap; th0c2(th0c2 >= theta_edges(end)) = theta_edges(end) - epsth;
ir_of_k = discretize(r0c2, r_edges);
it_of_k = discretize(th0c2, theta_edges);
if any(isnan(ir_of_k)) || any(isnan(it_of_k))
    error('k->(ir,it) mapping has NaN.');
end

% ===== 大网格材料：铜背景设为空气（空气带背景） =====
mat_code_big = mat_code;
mat_code_big(mat_code_big == 2) = 0;

% 大区并区键（空气/铁）
key_big = build_merge_keys_air_iron(mat_code_big);
keep_big_label = keep_one_label_per_region(nr, nt, key_big, ir_of_k, it_of_k);

% 铜岛并区键（空气/铁/铜分离），用于每个铜岛保留1个 label + 匝数
key_cu = build_merge_keys_air_iron_copper(mat_code);
keep_cu_label = keep_one_label_per_region(nr, nt, key_cu, ir_of_k, it_of_k);
comp_size = component_size_per_cell(nr, nt, key_cu, ir_of_k, it_of_k);

% 预构建铜布尔网格（基域）用于判断某个边是不是"铜岛外边界"
is_copper_grid = false(nr, nt);
for k = 1:Nd
    if mat_code(k) == 2
        is_copper_grid(ir_of_k(k), it_of_k(k)) = true;
    end
end

for s = 1:nSector
    pid = phase_id_sector{s};

    Cs = c0 + (s-1)*th_span;
    isMirror = mod(s,2)==0;

    if ~isMirror
        th_c = Cs + d;
    else
        th_c = Cs - d;
    end
    th_r = th_c;

    xc = r0c .* cosd(th_c); yc = r0c .* sind(th_c);
    xr = r0r .* cosd(th_r); yr = r0r .* sind(th_r);

    % 背景（空气/铁）删内部边界
    delete_internal_boundaries_sector(r_edges, theta_edges, c0, Cs, isMirror, nr, nt, key_big, ir_of_k, it_of_k);
    remove_isolated_nodes_sector(r_edges, theta_edges, c0, Cs, isMirror, nr, nt, key_big, ir_of_k, it_of_k);

    for k = 1:Nd
        code = mat_code(k);

        % ===== 大格子标签：铁保留铁，其余(含铜背景)为空气 =====
        if keep_big_label(k)
            if code == 1
                matName = mats.iron;
            else
                matName = mats.air;
            end
            mi_addblocklabel(xr(k), yr(k));
            mi_selectlabel(xr(k), yr(k));
            mi_setblockprop(matName, 1, 0, '', 0, groupIdRing, 0);
            mi_clearselected();
        end

        if code ~= 2
            continue;
        end

        % ===== 铜岛几何：仅对外边界边做 inset，内部共享边不 inset =====
        add_copper_cell_boundary_island_inset(k, s, isMirror, Cs, c0, ...
            nr, nt, r_edges, theta_edges, alpha, alpha, groupIdCuGeom, ...
            ir_of_k, it_of_k, is_copper_grid);

        % ===== 铜标签：每个铜连通岛仅保留一个 =====
        if keep_cu_label(k)
            phaseId = pid(k);
            if phaseId == 0
                circCu = ''; turns = 0;
            else
                circCu = circNames{phaseId};
                turns = comp_size(k) * turns_per_circ(phaseId);
            end
            mi_addblocklabel(xc(k), yc(k));
            mi_selectlabel(xc(k), yc(k));
            mi_setblockprop(mats.copper, 1, 0, circCu, 0, groupIdCore, turns);
            mi_clearselected();
        end
    end
end

end

% ======================================================================
% 给第 k 个铜格子画"岛级外边界内缩"小多边形
% 思路：某边若邻格也是铜 => 该边不内缩；否则按 alpha 内缩
% ======================================================================
function add_copper_cell_boundary_island_inset(k, s, isMirror, Cs, c0, ...
    nr, nt, r_edges, theta_edges, sh_r, sh_th, groupIdCuGeom, ...
    ir_of_k, it_of_k, is_copper_grid)

ir = ir_of_k(k);
it = it_of_k(k);

r1 = r_edges(ir);   r2 = r_edges(ir+1);
th1 = theta_edges(it); th2 = theta_edges(it+1);

dr = sh_r  * (r2 - r1);
dth = sh_th * (th2 - th1);

% 4邻域铜性（越界视为非铜）
has_in  = (ir > 1 ) && is_copper_grid(ir-1, it);
has_out = (ir < nr) && is_copper_grid(ir+1, it);
has_lft = (it > 1 ) && is_copper_grid(ir, it-1);
has_rgt = (it < nt) && is_copper_grid(ir, it+1);

% 外边界才内缩；内部共享边保持原位
ri1 = r1 + (~has_in ) * dr;
ri2 = r2 - (~has_out) * dr;
ti1 = th1 + (~has_lft) * dth;
ti2 = th2 - (~has_rgt) * dth;

% 防御：极端情况下避免反转
if ri2 <= ri1
    rm = 0.5*(r1+r2);
    ri1 = rm - 0.1*(r2-r1);
    ri2 = rm + 0.1*(r2-r1);
end
if ti2 <= ti1
    tm = 0.5*(th1+th2);
    ti1 = tm - 0.1*(th2-th1);
    ti2 = tm + 0.1*(th2-th1);
end

mapTh = @(th_base) (~isMirror) * (Cs + (th_base - c0)) + ...
                   ( isMirror) * (Cs - (th_base - c0));
T1 = mapTh(ti1);
T2 = mapTh(ti2);

p = zeros(4,2);
p(1,:) = [ri1*cosd(T1), ri1*sind(T1)];
p(2,:) = [ri2*cosd(T1), ri2*sind(T1)];
p(3,:) = [ri2*cosd(T2), ri2*sind(T2)];
p(4,:) = [ri1*cosd(T2), ri1*sind(T2)];

for i = 1:4
    mi_addnode(p(i,1), p(i,2));
    mi_selectnode(p(i,1), p(i,2));
    mi_setnodeprop('', groupIdCuGeom);
    mi_clearselected();
end

for i = 1:4
    j = i + 1; if j == 5, j = 1; end
    mi_addsegment(p(i,1), p(i,2), p(j,1), p(j,2));
    mx = 0.5*(p(i,1)+p(j,1)); my = 0.5*(p(i,2)+p(j,2));
    mi_selectsegment(mx, my);
    mi_setsegmentprop('', 0, 1, 0, groupIdCuGeom);
    mi_clearselected();
end

end

function thw = wrap_to_span(th, th_start, th_end)
span = th_end - th_start;
thw = mod(th - th_start, 360);
thw(thw >= span & thw < (360-span)) = thw(thw >= span & thw < (360-span)) - 360;
thw = thw + th_start;
end

function keyMerge = build_merge_keys_air_iron(mat_code)
Nd = numel(mat_code);
keyMerge = zeros(Nd,1);
for k = 1:Nd
    if mat_code(k) == 1
        keyMerge(k) = 200; % iron
    else
        keyMerge(k) = 100; % air
    end
end
end

function keyMerge = build_merge_keys_air_iron_copper(mat_code)
Nd = numel(mat_code);
keyMerge = zeros(Nd,1);
for k = 1:Nd
    if mat_code(k) == 0
        keyMerge(k) = 100;
    elseif mat_code(k) == 1
        keyMerge(k) = 200;
    else
        keyMerge(k) = 300;
    end
end
end

function comp_size = component_size_per_cell(nr, nt, key, ir_of_k, it_of_k)
keyGrid = nan(nr, nt); kGrid = nan(nr, nt);
for k = 1:numel(key)
    keyGrid(ir_of_k(k), it_of_k(k)) = key(k);
    kGrid(ir_of_k(k), it_of_k(k)) = k;
end
comp_size = ones(numel(key),1);
visited = false(nr, nt);
for i = 1:nr
    for j = 1:nt
        if visited(i,j) || isnan(keyGrid(i,j)), continue; end
        target = keyGrid(i,j);
        stack = [i, j]; visited(i,j) = true;
        comp = zeros(0,2);
        while ~isempty(stack)
            ci = stack(end,1); cj = stack(end,2); stack(end,:) = [];
            comp = [comp; ci, cj]; %#ok<AGROW>
            nb = [ci-1,cj; ci+1,cj; ci,cj-1; ci,cj+1];
            for n = 1:4
                ni = nb(n,1); nj = nb(n,2);
                if ni<1 || ni>nr || nj<1 || nj>nt || visited(ni,nj), continue; end
                if keyGrid(ni,nj) == target
                    visited(ni,nj) = true;
                    stack = [stack; ni, nj]; %#ok<AGROW>
                end
            end
        end
        sz = size(comp,1);
        for t = 1:sz
            kk = kGrid(comp(t,1), comp(t,2));
            if ~isnan(kk), comp_size(kk) = sz; end
        end
    end
end
end

function keep_label = keep_one_label_per_region(nr, nt, key, ir_of_k, it_of_k)
keyGrid = nan(nr, nt); kGrid = nan(nr, nt);
for k = 1:numel(key)
    keyGrid(ir_of_k(k), it_of_k(k)) = key(k);
    kGrid(ir_of_k(k), it_of_k(k)) = k;
end
keep_label = false(numel(key), 1);
visited = false(nr, nt);
for i = 1:nr
    for j = 1:nt
        if visited(i,j), continue; end
        if isnan(keyGrid(i,j)), visited(i,j) = true; continue; end
        target = keyGrid(i,j);
        stack = [i, j]; visited(i,j) = true;
        rep_k = kGrid(i,j);
        while ~isempty(stack)
            ci = stack(end,1); cj = stack(end,2); stack(end,:) = [];
            nb = [ci-1,cj; ci+1,cj; ci,cj-1; ci,cj+1];
            for n = 1:4
                ni = nb(n,1); nj = nb(n,2);
                if ni<1 || ni>nr || nj<1 || nj>nt || visited(ni,nj), continue; end
                if keyGrid(ni,nj) == target
                    visited(ni,nj) = true;
                    stack = [stack; ni, nj]; %#ok<AGROW>
                end
            end
        end
        if ~isnan(rep_k), keep_label(rep_k) = true; end
    end
end
end

function delete_internal_boundaries_sector(r_edges, theta_edges, c0, Cs, isMirror, nr, nt, key, ir_of_k, it_of_k)
keyGrid = nan(nr, nt);
for k = 1:numel(key)
    keyGrid(ir_of_k(k), it_of_k(k)) = key(k);
end
for j = 1:nt
    th_mid_base = 0.5*(theta_edges(j) + theta_edges(j+1));
    th_mid = map_theta(th_mid_base, c0, Cs, isMirror);
    for i = 1:(nr-1)
        if keyGrid(i,j) == keyGrid(i+1,j)
            r = r_edges(i+1); x = r*cosd(th_mid); y = r*sind(th_mid);
            try
                mi_clearselected(); mi_selectarcsegment(x,y); mi_deleteselected(); mi_clearselected();
            catch
                try
                    mi_clearselected(); mi_selectsegment(x,y); mi_deleteselected(); mi_clearselected();
                catch
                    mi_clearselected();
                end
            end
        end
    end
end
for j = 1:(nt-1)
    th_base = theta_edges(j+1);
    th = map_theta(th_base, c0, Cs, isMirror);
    for i = 1:nr
        if keyGrid(i,j) == keyGrid(i,j+1)
            rmid = 0.5*(r_edges(i) + r_edges(i+1));
            x = rmid*cosd(th); y = rmid*sind(th);
            try
                mi_clearselected(); mi_selectsegment(x,y); mi_deleteselected(); mi_clearselected();
            catch
                mi_clearselected();
            end
        end
    end
end
end

function remove_isolated_nodes_sector(r_edges, theta_edges, c0, Cs, isMirror, nr, nt, key, ir_of_k, it_of_k)
keyGrid = nan(nr, nt);
for k = 1:numel(key)
    keyGrid(ir_of_k(k), it_of_k(k)) = key(k);
end
for i = 2:nr
    r = r_edges(i);
    for j = 2:nt
        k11 = keyGrid(i-1, j-1); k21 = keyGrid(i,   j-1);
        k12 = keyGrid(i-1, j);   k22 = keyGrid(i,   j);
        if isnan(k11) || isnan(k21) || isnan(k12) || isnan(k22), continue; end
        if ~(k11 == k21 && k11 == k12 && k11 == k22), continue; end
        th_base = theta_edges(j);
        th = map_theta(th_base, c0, Cs, isMirror);
        x = r*cosd(th); y = r*sind(th);
        try
            mi_clearselected(); mi_selectnode(x,y); mi_deleteselected(); mi_clearselected();
        catch
            mi_clearselected();
        end
    end
end
end

function th_out = map_theta(th_in, c0, Cs, isMirror)
d = th_in - c0;
if ~isMirror
    th_out = Cs + d;
else
    th_out = Cs - d;
end
end
