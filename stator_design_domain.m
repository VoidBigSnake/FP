function domain = stator_design_domain(cfg)
%STATORDESIGNDOMAIN Build a 15-degree stator-sector design mask.
%  DOMAIN = STATORDESIGNDOMAIN(CFG) returns grid edges and logical masks
%  describing which polar cells can be edited by the optimizer.
%
%  Required fields in CFG:
%    nr             - number of radial cells (e.g., 8)
%    nt             - number of tangential cells within 15° (e.g., 20)
%    r_inner        - stator inner radius (airgap side)
%    r_outer        - stator yoke outer radius
%    slot_r_inner   - slot window inner radius (start of slot opening)
%    slot_r_outer   - slot window outer radius (end of slot opening)
%    slot_span_deg  - angular span of the slot window within the 15° sector
%    theta_span_deg - total angular size of the design sector (15)
%
%  Optional fields:
%    coil_keepout_deg - extra half-angle to freeze around the slot mouth
%                       (default 0)
%    yoke_buffer_deg  - half-angle buffer near sector edges to keep fixed
%                       for better periodic boundaries (default 0)
%
%  Output DOMAIN fields:
%    r_edges, theta_edges - grid edges
%    design_mask          - nt-by-nr logical, true where optimization is
%                           allowed (iron/void toggle); false means frozen
%    slot_mask            - nt-by-nr logical, true for the slot/coil region
%                           that must stay fixed
%    notes                - cell array with textual reminders
%
%  The intent is to keep the slot window and coil region immutable while
%  letting the optimizer edit iron in the yoke and tooth-shoulder area.
%
%  Example:
%    cfg = struct('nr',8,'nt',20,'r_inner',22,'r_outer',40,...
%                 'slot_r_inner',22.5,'slot_r_outer',28,...
%                 'slot_span_deg',6,'theta_span_deg',15,...
%                 'coil_keepout_deg',1.0,'yoke_buffer_deg',0.5);
%    domain = stator_design_domain(cfg);
%
%  See also APPLY_DESIGN_MASK.

arguments
    cfg.nr (1,1) double {mustBePositive}
    cfg.nt (1,1) double {mustBePositive}
    cfg.r_inner (1,1) double {mustBePositive}
    cfg.r_outer (1,1) double {mustBePositive}
    cfg.slot_r_inner (1,1) double {mustBePositive}
    cfg.slot_r_outer (1,1) double {mustBePositive}
    cfg.slot_span_deg (1,1) double {mustBePositive}
    cfg.theta_span_deg (1,1) double {mustBePositive}
    cfg.coil_keepout_deg (1,1) double {mustBeNonnegative} = 0
    cfg.yoke_buffer_deg (1,1) double {mustBeNonnegative} = 0
end

if cfg.r_outer <= cfg.r_inner
    error('r_outer must be greater than r_inner (got %.2f vs %.2f).', ...
          cfg.r_outer, cfg.r_inner);
end

if cfg.slot_r_outer <= cfg.slot_r_inner
    error('slot_r_outer must be greater than slot_r_inner (got %.2f vs %.2f).', ...
          cfg.slot_r_outer, cfg.slot_r_inner);
end

% Build edges and centers.
r_edges = linspace(cfg.r_inner, cfg.r_outer, cfg.nr + 1);
theta_edges = linspace(0, cfg.theta_span_deg, cfg.nt + 1); % degrees
[theta_centers, r_centers] = meshgrid(
    theta_edges(1:end-1) + diff(theta_edges)/2,
    r_edges(1:end-1) + diff(r_edges)/2);
% theta_centers is nt-by-1 replicated across rows; transpose for masks.
theta_centers = theta_centers';
r_centers = r_centers';

% Slot window mask (tangential span around the center of the sector).
slot_half = cfg.slot_span_deg/2 + cfg.coil_keepout_deg;
slot_center = cfg.theta_span_deg/2;
theta_slot = abs(theta_centers - slot_center) <= slot_half;
r_slot = (r_centers >= cfg.slot_r_inner) & (r_centers <= cfg.slot_r_outer);
slot_mask = theta_slot & r_slot;

% Buffer near sector edges to keep periodic boundary nodes intact.
edge_mask = (theta_centers <= cfg.yoke_buffer_deg) | ...
            (theta_centers >= (cfg.theta_span_deg - cfg.yoke_buffer_deg));

% Design mask: allow editing only in yoke/tooth shoulder, not slots/edges.
design_mask = ~slot_mask & ~edge_mask;

% Package result.
domain = struct();
domain.r_edges = r_edges;
domain.theta_edges = theta_edges;
domain.design_mask = design_mask;
domain.slot_mask = slot_mask;
domain.notes = {
    'design_mask: true = optimizer can toggle iron/void';
    'slot_mask: true = slot/coil region kept immutable';
    'apply_design_mask(bits, domain, base_val) enforces frozen regions'
    sprintf('sector: %.1f deg, slot span: %.1f deg (keepout %.1f deg)', ...
            cfg.theta_span_deg, cfg.slot_span_deg, cfg.coil_keepout_deg)
};

end

function grid = apply_design_mask(bits, domain, base_val)
%APPLY_DESIGN_MASK Reshape genome bits onto the grid and freeze slot cells.
%  GRID = APPLY_DESIGN_MASK(BITS, DOMAIN, BASE_VAL) returns an nt-by-nr
%  array where editable cells come from BITS (reshaped), and frozen cells
%  (slot/coil or edge buffers) are forced to BASE_VAL (e.g., iron = 1).
%
%  BITS is expected to be a vector of length nt*nr matching DOMAIN masks.

if numel(bits) ~= numel(domain.design_mask)
    error('Bitstring length %d does not match grid size %d.', ...
          numel(bits), numel(domain.design_mask));
end

grid = reshape(bits, size(domain.design_mask));
freeze_mask = ~domain.design_mask; % includes slot and edge buffers
grid(freeze_mask) = base_val;
end
