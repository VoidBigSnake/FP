function [J, T_avg, T_ripple] = eval_design_femm_worker(bits, ctx, S)

    domain          = ctx.domain;
    phase_id_sector = ctx.phase_id_sector;
    N_phase_total   = ctx.N_phase_total;

    [mat_code, turns_per_cell] = compute_turns_per_cell(bits, domain, phase_id_sector, N_phase_total);

    if all(turns_per_cell == 0)
        J = 1e6; T_avg = 0; T_ripple = 0;
        return;
    end

    % 关键：不要 opendocument 每次都来一遍
    % 只需要删旧设计块，然后重新打 label
    mi_selectgroup(ctx.groupId_core);
    mi_deleteselected();
    mi_selectgroup(ctx.groupId_ring);
    mi_deleteselected();

    femm_apply_design_bits_rep6_inset(bits, domain,ctx, ctx.phase_id_sector, ctx.mats, ctx.circNames, ...
                             ctx.groupId_core, ctx.groupId_ring,31, turns_per_cell);
  
    % 强烈建议：每次分析前把模型保存到 worker 私有目录、并且文件名唯一
% --- 每次评估保存到唯一case（并行安全） ---
femCase = fullfile(S.tmpDir, ['case_' char(java.util.UUID.randomUUID) '.fem']);
mi_saveas(femCase);

    % scantorque 里如果自己也 mi_saveas / 导出文件，请同样改成用 femCase 前缀
    [T_avg, T_ripple] = scantorque(3.5);

    J_ripple = T_ripple;
    if T_avg < ctx.T_min
        penalty = ctx.penaltyCoef * (ctx.T_min - T_avg)/ctx.T_min;
    else
        penalty = 0;
    end
    J = J_ripple + penalty;
end
