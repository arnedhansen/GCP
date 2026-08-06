function [info_lines, info_viol] = build_selection_gate_info_columns( ...
    ci, occdom_vec, emg_temp_vec, emg_hf_slope_vec)
% Build QC text lines and fail flags for occdom, EMG_temp, and EMG_hf_slope.
occdom_val = occdom_vec(ci);
emg_temp_val = emg_temp_vec(ci);
emg_hf_val = emg_hf_slope_vec(ci);
fail_occdom = ~(isfinite(occdom_val) && occdom_val > 1);
fail_emg_temp = isfinite(emg_temp_val) && emg_temp_val >= 1;
fail_emg_hf = isfinite(emg_hf_val) && emg_hf_val > 0;
info_lines = { ...
    sprintf('occdom: %.2f (> 1)', occdom_val), ...
    sprintf('EMG_temp: %.2f (< 1)', emg_temp_val), ...
    sprintf('EMG_hf_slope: %.2f (<= 0)', emg_hf_val)};
info_viol = [fail_occdom, fail_emg_temp, fail_emg_hf];
end
