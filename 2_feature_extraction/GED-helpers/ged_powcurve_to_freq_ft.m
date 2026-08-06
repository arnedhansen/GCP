function freq = ged_powcurve_to_freq_ft(pow_vec, scan_freqs, data_recording)
% Wrap a GED power curve as a one-channel FieldTrip freq structure.
% Does not set topolabel (FieldTrip would treat the struct as a component).
freq = struct();
freq.label = {'GED'};
freq.freq = scan_freqs(:)';
freq.dimord = 'chan_freq';
nf = numel(scan_freqs);
if isempty(scan_freqs)
    freq.powspctrm = nan(1, 0);
    return;
end
if isempty(pow_vec) || ~isnumeric(pow_vec) || numel(pow_vec) ~= nf
    freq.powspctrm = nan(1, nf);
else
    freq.powspctrm = reshape(double(pow_vec(:)).', 1, nf);
end
if nargin >= 3 && isstruct(data_recording)
    if isfield(data_recording, 'elec')
        el = data_recording.elec;
        if isstruct(el) && isfield(el, 'label') && ~isempty(el.label)
            freq.elec = el;
        end
    end
    if isfield(data_recording, 'grad')
        g = data_recording.grad;
        if isstruct(g) && isfield(g, 'label') && ~isempty(g.label)
            freq.grad = g;
        end
    end
end
end
