### README for GCP Study (Gamma Contrast Perception)

Combined EEG and eye-tracking (ET) study of gamma oscillations and oculomotor dynamics during perception of dynamic concentric inward gratings at four contrast levels (25%, 50%, 75%, 100%).

The titles below correspond to the folder names. Apart from the Python script in `4_stats`, all files are MATLAB scripts.

## paradigms

The grating task is run with `master.m`, which calls `GCP_gratingsTask.m` (PsychToolbox). Each session includes resting EEG, a short training block, and four experimental blocks (704 trials total: 176 per contrast). Triggers 61–64 code the four contrast conditions (no-task blocks); triggers 51–54 code task blocks with button responses. Dependencies (EEG/ET setup, calibration, luminance conversion, screen settings) are in the `paradigm` folder.

## 1_preprocessing

### 1_cut

Raw ANT `.cnt` EEG files are cut into resting and block recordings by `GCP_doCutting.m` using `GCP_cutData.m`. EyeLink `.asc` files are converted and synchronized with the cut EEG in the same script.

### 2_automagic

Cut EEG and ET files are preprocessed in Automagic (Pedroni et al., 2019). Bad channels are detected with the EEGLAB plugin clean_rawdata (Mullen et al., 2015). Data are high-pass filtered at 0.1 Hz, line noise is removed with ZapLine (De Cheveigné, 2020), and ocular artifacts are corrected with OPTICAT (Dimigen, 2020) and ICLabel (Pion-Tonachini et al., 2019). Bad electrodes are interpolated. Blocks failing Automagic quality criteria are excluded.

### 3_merge

Automagic-preprocessed EEG is merged with the corresponding ET files using `GCP_mergeData.m` (EYE-EEG `pop_importeyetracker`).

### 4_preprocessing

`GCP_preprocessing.m` segments merged EEG+ET data into epochs around stimulus onset ([-2, 3.5] s), converts to FieldTrip format, and writes per-condition EEG and ET files. Baseline for gaze metric extraction: [-1.5, -0.5] s; analysis window for scalar gaze metrics: [0.3, 2] s. EyeLink event counts (blinks, fixations, saccades) are extracted per trial during preprocessing. Outputs go to `data/features/<subject>/eeg` (`dataEEG.mat`) and `.../gaze` (`dataET.mat`, `gaze_metrics.mat`).

## 2_feature_extraction

`GCP_behavioral_fex.m` extracts accuracy and reaction time per trial and per subject-condition mean. Outputs: `behavioral_matrix_trial.mat`, `behavioral_matrix_subj.mat`, and group-level `GCP_behavioral_matrix.mat`.

`GCP_gaze_fex.m` extracts gaze deviation, gaze SD, BCEA (k = 2.291, 95%), pupil size, microsaccade rate, and eye velocity (horizontal, vertical, 2D). Scalar metrics are computed for baseline [-1.5, -0.5] s and stimulus [0, 2] s, with dB baseline correction (`10*log10(stim/baseline)`). Time courses for pupil, velocity, and microsaccades are saved for visualization. Outputs: `gaze_matrix_trial.mat`, `gaze_matrix_subj.mat`, and group-level `GCP_gaze_matrix.mat`.

`GCP_eeg_fex_GED.m` is the primary EEG analysis. Per subject, trials are pooled across conditions and gamma-band (30–90 Hz) covariances are computed for baseline [-1.5, -0.5] s and three stimulus windows (full [0, 2] s, early [0, 0.5] s, late [1, 2] s). Window-specific generalized eigendecomposition (GED) is solved with regularization; candidate components are ranked by eigenvalue and scored on occipital topography, spectral form, and artifact metrics. An eigenvalue-weighted combined component is built per subject. Each trial is projected to this component space and scanned on a 30–90 Hz grid (mtmfft, 3 Hz multitaper smoothing). Per-trial peak gamma frequency and peak power (mean power within peak ± 5 Hz) are extracted; unstable trials are flagged automatically. Outputs: `GCP_eeg_GED.mat` (trial-level and subject-level metrics, component diagnostics, outlier masks) and `GCP_eeg_powspctrm_GED.mat` (FieldTrip freq structs for grand-average spectra). Subject inclusion for downstream GED analyses is written to `controls/GCP_subject_inclusion.mat` (subjects with valid gamma power).

`GCP_eeg_fex_GED_TFR.m` reconstructs subject-specific GED spatial filters, projects single-trial EEG, and computes baseline-corrected TFRs (30–90 Hz, ERS/ERD-style dB subtraction). Output: `GCP_eeg_GED_TFR.mat`.

`GCP_master_matrix.m` merges subject-level behavioral, gaze, and GED tables into `GCP_merged_data.mat` / `.csv`. GED columns use median peak frequency and robust mean peak power from `GCP_eeg_GED.mat`.

`GCP_master_matrix_trials.m` merges trial-level behavioral, gaze, and GED tables into `GCP_merged_data_trials.mat` / `.csv`. Both master matrices attach the `Include` flag from `GCP_subject_inclusion.mat`.

**Run order:** `4_preprocessing` → behavioral → gaze → GED → GED TFR → master matrices. The trial-level CSV is the input for the Python raincloud script.

## 3_visualization

**Behavioral:** Accuracy by condition (`behavioral/GCP_behav.m`).

**EEG (GED):** Grand-average and single-subject GED power spectra (`eeg/powspctrm/GCP_eeg_powspctrm_GED.m`); GED-projected TFRs and 100% − 25% difference maps (`eeg/tfr/GCP_TFR_GED.m`); pooled trial-level boxplots for gamma peak frequency and peak power (`eeg/GCP_eeg_GED_trial_boxplots.m`).

**Gaze:** Baseline-normalized time courses with SEM shading for microsaccades (`gaze/microsaccades/GCP_gaze_microsaccades_TC.m`), pupil size (`gaze/pupilSize/GCP_gaze_pupilSize_TC.m`), eye velocity (`gaze/velocity/GCP_gaze_velocity_TC.m`), and fixation rate reconstructed from EyeLink fixation onsets (`gaze/fixations/GCP_gaze_fixations_TC.m`). Microsaccade, pupil, and velocity plots read dB-baselined time courses saved by `GCP_gaze_fex.m`.

**Hypotheses:** Schematic gamma spectra by contrast (`hypotheses/GCP_hypotheses_plot.m`).

All visualization scripts read from `data/features/` and write figures to `figures/`. GED-related plots apply the subject inclusion list via `gcp_subject_inclusion`. Run after feature extraction.

## 4_stats

### Subject-level overview (MATLAB)

`GCP_stats_overview.m` loads `GCP_merged_data.mat` and produces a multi-panel overview of all numeric variables by contrast condition. `GCP_stats_boxplots.m` produces one figure per variable with boxplots, subject lines, and jittered dots (GED cohort filter applied when `Include` is present).

### Hypothesis testing (MATLAB)

`GCP_hypotheses_trials.m` tests registered hypotheses on trial-level data with raincloud-style plots:

- **Oculomotor (gaze):** H1 microsaccade rate decreases with contrast; H2 eye velocity increases with contrast; H3 pupil constriction amplitude increases with contrast; H4 BCEA increases with contrast.
- **Gamma (GED):** H5 gamma peak frequency increases with contrast; H6 gamma peak amplitude peaks at 75% contrast (inverted U).
- **Cross-modal:** H7 gamma frequency relates to oculomotor dynamics (subject-level means, since gaze and EEG trials are not matched trial-by-trial).

Trial outliers are excluded data-driven (median ± 3 MAD per variable).

### Rainclouds (Python)

`GCP_stats_rainclouds.py` produces trial-level raincloud figures (half-kernel densities, boxplots, jittered trials) for gamma peak frequency, gamma peak power, microsaccade rate (dB), pupil size (dB), gaze velocity Y (dB), and reaction time. GED trial metrics are reconstructed from `GCP_eeg_GED.mat` (same peak-power definition as the MATLAB pipeline). Inputs: `GCP_merged_data_trials.csv` and `GCP_eeg_GED.mat`. Optional MixedLM significance brackets are available but off by default. Python helpers (`stats_helpers`, `rainclouds_plotting_helpers`) come from [github.com/arnedhansen/functions](https://github.com/arnedhansen/functions). Adapt `base_dir` and input paths in the script to your setup.

## Additional Files

### controls

`GCP_subject_inclusion.m` defines a manual seed/override for the analysis cohort (default exclusion of subjects 602, 604, 608). After `GCP_eeg_fex_GED.m` runs, inclusion is overwritten automatically based on valid gamma power. `GCP_et_cal_val.m` summarises EyeLink calibration and validation quality. `GCP_eeg_impedances.m` extracts and plots ANT impedances before and after recording.

### tests

Exploratory GED–microsaccade coupling analyses (`GCP_test_GED_microsaccade_*.m`) and replication scripts (`GCP_Replication_Kuo_Juan.m`).

### Dependencies

`startup` and `setup('GCP')` (paths, subject list, colours, head model; `setup` is in [github.com/arnedhansen/functions](https://github.com/arnedhansen/functions)). FieldTrip and EEGLAB are required for preprocessing and spectral analysis. For plots: `shadedErrorBar` (time courses), `color_def('GCP')` (condition colours). Many scripts hardcode data roots (e.g. `/Volumes/g_psyplafor_methlab$/Students/Arne/GCP` or `W:\...`); change these to your `data/` location.
