### README for GCP Study (Gamma Contrast Perception)

Combined EEG and eye-tracking (ET) study of gamma oscillations and oculomotor dynamics during perception of dynamic concentric inward gratings at four contrast levels (25%, 50%, 75%, 100%).

The titles below correspond to the folder names. Apart from the Python scripts in `4_stats`, all files are MATLAB scripts.

## paradigms

The grating task is run with `master.m`, which calls `GCP_gratingsTask.m` (PsychToolbox). Each session includes resting EEG, a short training block, and four experimental blocks (704 trials total: 176 per contrast). Triggers 61 to 64 code the four contrast conditions (no-task trials). Triggers 51 to 54 code catch trials with a white fixation and a space-bar response. Dependencies (EEG/ET setup, calibration, luminance conversion, screen settings) are in the `paradigm` folder.

## 1_preprocessing

### 1_cut

Raw ANT `.cnt` EEG files are cut into resting and block recordings by `GCP_doCutting.m` using `GCP_cutData.m`. EyeLink `.asc` files are converted and synchronized with the cut EEG in the same script.

### 2_automagic

Cut EEG and ET files are preprocessed in Automagic (Pedroni et al., 2019). Bad channels are detected with the EEGLAB plugin clean_rawdata (Mullen et al., 2015). Data are high-pass filtered at 0.1 Hz, line noise is removed with ZapLine (De Cheveigné, 2020), and ocular artifacts are corrected with OPTICAT (Dimigen, 2020) and ICLabel (Pion-Tonachini et al., 2019). Bad electrodes are interpolated. Blocks failing Automagic quality criteria are excluded.

### 3_merge

Automagic-preprocessed EEG is merged with the corresponding ET files using `GCP_mergeData.m` (EYE-EEG `pop_importeyetracker`).

### 4_preprocessing

`GCP_preprocessing.m` opens the AOC `askRunMode` dialog (ALL subjects or only NEW subjects without `dataEEG.mat`; Cancel aborts). It segments merged EEG+ET data into epochs around stimulus onset ([−2, 3.5] s), converts to FieldTrip, and writes per-condition EEG and ET files.

Among non-catch trials (triggers 61 to 64), failed pre-stimulus fixation (`fixation == 0`) is dropped. Replacements with `fixation == 1` are kept. Epoch count is asserted against the block behavioural mat.

EEG is average-rereferenced, then epochs with any channel beyond ±90 μV are rejected (`ft_artifact_threshold` / `ft_rejectartifact`). If more than 50% of the post-fixation epochs in a block are flagged, the whole block is discarded.

Baseline for gaze metric extraction: [−1.5, −0.5] s. Confirmatory analysis window for scalar gaze metrics: [0, 2] s. EyeLink blink, fixation, and saccade rates are counted on the same remaining epochs as saved `dataET`. Outputs go to `data/features/<subject>/eeg` (`dataEEG.mat`) and `.../gaze` (`dataET.mat`, `gaze_metrics.mat`).

## 2_feature_extraction

`GCP_behavioral_fex.m` extracts catch-hit accuracy (space on white-fixation trials / all white-fixation trials) and reaction time. `WhiteCross` is kept on the trial table. Outputs: `behavioral_matrix_trial.mat`, `behavioral_matrix_subj.mat`, and group-level `GCP_behavioral_matrix.mat`.

`GCP_gaze_fex.m` extracts BCEA (k = 5.991, 95%), pupil size, microsaccade rate, eye velocity, and EyeLink event counts. Blinks are invalid/missing gaze and/or pupil at floor, padded ±50 ms (25 samples at 500 Hz). Microsaccades are Engbert events with 2D amplitude 0.1 to 1.0° at analysis `ppd = 50`. Baselined (% change) scalars use suffix `_bl` (`MSRate_bl`, `BCEA_bl`, …) for the full [0, 2] s window. Subject×condition confirmatory scalars are the mean of trial-level percentage change for microsaccades, BCEA, pupil, and velocity. Group file: `GCP_gaze_window_summaries.mat`.

`GCP_eeg_fex_GED.m` is the primary EEG analysis. Per subject, trials are pooled across contrast and gamma-band (30 to 90 Hz) covariances are computed for baseline [−1.5, −0.5] s and the full stimulus window [0, 2] s. Generalized eigendecomposition (GED) is solved with regularization; candidate components are ranked by eigenvalue and scored on occipital topography, spectral form, and artifact metrics. An eigenvalue-weighted combined component is built per subject. Each trial is projected to this component space and scanned on a 30 to 90 Hz grid (mtmfft, 3 Hz multitaper smoothing). Per-trial peak gamma frequency and peak power (mean power within peak ± 5 Hz) are extracted; unstable trials are flagged automatically. Outputs: `GCP_eeg_GED.mat` and `GCP_eeg_powspctrm_GED.mat`. Subject inclusion is written to `controls/GCP_subject_inclusion.mat`. Optional GED-projected TFRs go to `GCP_eeg_GED_TFR.mat`.

`GCP_master_matrix.m` merges subject-level behavioral, gaze, and GED tables into `GCP_merged_data.mat` / `.csv` (Include flag from `GCP_subject_inclusion.mat`). Confirmatory GED columns are full-window `Power` and `Frequency`. Confirmatory gaze columns include `PupilSize_bl`, `MSRate_bl`, `BCEA_bl`, and `Vel2D_bl`.

**Run order:** `4_preprocessing` → behavioral → gaze → GED → master matrix.

## 3_visualization

**Behavioral:** Accuracy by condition (`behavioral/GCP_behav.m`).

**EEG (GED):** Grand-average and single-subject GED power spectra (`eeg/powspctr/GCP_eeg_powspctrm_GED.m`); GED-projected TFRs with pairwise CBPT cluster outlines (`eeg/tfr/GCP_TFR_GED.m`; 2×2 GA, 100% vs 25%, and all six pairwise difference maps); pooled trial-level boxplots for gamma peak frequency and peak power (`eeg/GCP_eeg_GED_trial_boxplots.m`).

**Gaze:** Baseline-normalized time courses with SEM shading for microsaccades (`gaze/microsaccades/GCP_gaze_microsaccades_TC.m`), pupil size (`gaze/pupilSize/GCP_gaze_pupilSize_TC.m`), eye velocity (`gaze/velocity/GCP_gaze_velocity_TC.m`), and fixation rate reconstructed from EyeLink fixation onsets (`gaze/fixations/GCP_gaze_fixations_TC.m`). Microsaccade, pupil, and velocity plots read dB-baselined time courses saved by `GCP_gaze_fex.m`. Full-window BCEA ellipses: `gaze/bcea/GCP_gaze_BCEA_ellipses.m`.

**Hypotheses:** Schematic gamma spectra by contrast (`hypotheses/GCP_hypotheses_plot.m`).

All visualization scripts read from `data/features/` and write figures to `figures/`. GED-related plots apply the subject inclusion list via `gcp_subject_inclusion`. Run after feature extraction.

## 4_stats

### Subject-level overview (MATLAB)

`GCP_stats_overview.m` loads `GCP_merged_data.mat` and produces a multi-panel overview of all numeric variables by contrast condition. `GCP_stats_boxplots.m` loads precomputed full-window subject×condition scalars from `GCP_eeg_GED.mat` and `GCP_gaze_window_summaries.mat`. Output: `figures/stats/boxplots/`.

### Confirmatory LMMs (Python)

`GCP_stats_lmm.py` fits full-window subject×contrast MixedLM models on `GCP_merged_data.csv` (`Include == 1`):

- **Primary (H1 to H6):** continuous linear contrast slope (`contrast_num_c`) for pupil, microsaccade rate, BCEA, eye velocity (Vel2D), gamma peak power, and gamma peak frequency.
- **Follow-up:** categorical MixedLM with all six pairwise contrasts, FDR-BH within each DV.
- **H7:** `gamma ~ gaze_c * contrast_num_c` for Frequency and Power each paired with `MSRate_bl`, `BCEA_bl`, and `Vel2D_bl`. `gaze_c` is centred within subject. Confirmatory term is the main effect of `gaze_c`. FDR-BH across the six gaze main-effect p-values (`GCP_mixedlm_h7.csv`).

Writes AOC-style CSVs under `data/stats/` (`GCP_mixedlm_slope.csv`, `GCP_pairwise_mixedlm.csv`, `GCP_mixedlm_h7.csv`, per-DV fixed tables).

### Rainclouds (Python)

`GCP_stats_rainclouds.py` plots subject-level full-window rainclouds for the same confirmatory DVs. FDR pairwise asterisks annotate the figures; the primary slope is printed to the console. Prefers pairwise p-values from `GCP_stats_lmm.py` output when available. Python helpers (`stats_helpers`, `rainclouds_plotting_helpers`) come from [github.com/arnedhansen/functions](https://github.com/arnedhansen/functions).

### Secondary CBPT (MATLAB)

`GCP_stats_ftests.m`: omnibus FieldTrip CBPT over time for pupil, microsaccade rate, and eye velocity. `GCP_TFR_GED.m`: pairwise GED TFR CBPT. These are secondary, not the confirmatory slope tests.

## Additional Files

### controls

`GCP_subject_inclusion.m` defines a manual seed/override for the analysis cohort (default exclusion of subjects 602, 604, 608). After `GCP_eeg_fex_GED.m` runs, inclusion is overwritten automatically based on valid gamma power. `GCP_et_cal_val.m` summarises EyeLink calibration and validation quality. `GCP_eeg_impedances.m` extracts and plots ANT impedances before and after recording.

### tests

Exploratory GED–microsaccade coupling analyses (`GCP_test_GED_microsaccade_*.m`) and replication scripts (`GCP_Replication_Kuo_Juan.m`).

### Dependencies

`startup` and `setup('GCP')` (paths, subject list, colours, head model; `setup` is in [github.com/arnedhansen/functions](https://github.com/arnedhansen/functions)). FieldTrip and EEGLAB are required for preprocessing and spectral analysis. For plots: `shadedErrorBar` (time courses), `color_def('GCP')` (condition colours). Many scripts hardcode data roots (e.g. `/Volumes/g_psyplafor_methlab$/Students/Arne/GCP` or `W:\...`); change these to your `data/` location.
