function [info_lines, info_viol] = build_selection_gate_info_columns( ...
    ci, post_front_vec, post_temp_vec, topo_peak_frac_vec, emg_hf_slope_vec)
% Build QC text lines and fail flags for spatial/EMG selection gates.
post_front_val = post_front_vec(ci);
post_temp_val = post_temp_vec(ci);
topo_peak_frac_val = topo_peak_frac_vec(ci);
emg_hf_val = emg_hf_slope_vec(ci);
fail_post_front = ~(isfinite(post_front_val) && post_front_val > 1);
fail_post_temp = ~(isfinite(post_temp_val) && post_temp_val > 1);
fail_topo_peak = isfinite(topo_peak_frac_val) && topo_peak_frac_val > 0.25;
fail_emg_hf = isfinite(emg_hf_val) && emg_hf_val > 0;
info_lines = { ...
    sprintf('post_front: %.2f (> 1)', post_front_val), ...
    sprintf('post_temp: %.2f (> 1)', post_temp_val), ...
    sprintf('topo_peak_frac: %.2f (<= 0.25)', topo_peak_frac_val), ...
    sprintf('EMG_hf_slope: %.2f (<= 0)', emg_hf_val)};
info_viol = [fail_post_front, fail_post_temp, fail_topo_peak, fail_emg_hf];
end
