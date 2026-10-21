%% GCP microsaccade locked low frequency ITPC (placeholder)
%
% Exploratory follow up. Not confirmatory. Do not run until implemented.
% Companion to GCP_ms_locked_GED_gamma.m, which fits a peri microsaccade
% spatial filter against a no microsaccade window. This script is the
% phase test from Leszczynski and Schroeder (2019) and is still a placeholder.
%
% Question
%   Does theta / alpha phase (4 to 12 Hz) cluster at Engbert microsaccade
%   onset, more than at stimulus onset or at jittered surrogate times?
%
% Rationale
%   The paper's mechanism is predictive reset of ongoing low frequency
%   oscillations to a high excitability phase at movement onset. Gamma
%   power (script 1) can change without that reset. Inter trial phase
%   clustering (ITPC) at MS onset is the direct test of the mechanism.
%   Stronger ITPC locked to saccade / microsaccade onset than to fixation
%   or stimulus onset is what they treat as a motor, not retinal, signal.
%
% Prediction
%   Occipital 4 to 12 Hz ITPC peaks at microsaccade onset and exceeds a
%   timing surrogate. Stimulus locked ITPC at grating onset is a separate
%   evoked alignment and is not evidence for the motor reset claim.
%
% Inputs (already on disk after preprocessing and gaze fex)
%   paths.features/<subj>/eeg/dataEEG.mat
%       dataEEG_c25 / c50 / c75 / c100
%   paths.features/<subj>/gaze/gaze_microsaccade_events.mat
%       ms_events_*.Onset{trl} in seconds on the trial time axis
%   controls/GCP_subject_inclusion.mat via gcp_subject_inclusion
%
% Spatial filter
%   Do not use the gamma GED filter from GCP_eeg_GED.mat. That filter was
%   fit on 30 to 90 Hz covariance and is the wrong spatial weights for
%   theta / alpha phase. Use occipital channels (O / PO / I) via
%   ft_selectdata, or a Laplacian on those sites (ft_scalpcurrentdensity)
%   if volume conduction of the evoked field is a concern.
%
% FieldTrip (do not wrap these in custom functions)
%   ft_selectdata              occipital channel subset
%   ft_redefinetrial           t = 0 at each microsaccade onset
%   ft_freqanalysis            output fourier (wavelet or mtmconvol)
%                              keeptrials = yes, foi 4 to 12 Hz
%   ITPC from the trial Fourier coefficients:
%       itpc = abs(mean(exp(1i * angle(freq.fourierspctrm)), 1))
%   Repeat the same ft_freqanalysis on stimulus locked epochs and on
%   jittered MS times. No custom phase toolbox required.
%
% Pipeline (implement later)
%   1. Load inclusion cohort, then per subject EEG and MS onsets.
%   2. Select occipital channels. Do not project through gamma GED.
%   3. Keep onsets in 0 to 2 s with enough padding for the wavelet /
%      taper at 4 Hz (low frequencies need long windows).
%   4. ft_redefinetrial: one epoch per retained MS, time 0 = onset.
%   5. Fourier time frequency, keep trials. Compute ITPC across MS epochs
%      within subject, then average ITPC maps across subjects.
%   6. Controls, same subject, same epoch length:
%        a) stimulus locked (t = 0 at grating onset)
%        b) surrogate (MS times jittered or circularly shifted)
%      Motor reset predicts MS locked ITPC > surrogate. Stimulus locked
%      ITPC can be high for a different reason (evoked response).
%   7. Optional: ITPC at MS onset vs a later sample (putative fixation /
%      reafference). The paper argues reset is stronger at movement onset.
%   8. Contrast split is secondary. Pooling events is safer because high
%      contrast yields fewer microsaccades.
%
% Caveats
%   ITPC rises with epoch count. Report N events per subject and match
%   surrogate N to real N. Do not compare conditions with unequal event
%      counts without that control.
%   OPTICAT / ICLabel can remove genuine peri movement cortical phase.
%   If implemented, test whether occipital ITPC at MS onset depends on
%   that ICA step.
%   Instructed fixation plus a 2 Hz drift is not exploratory saccadic
%   sampling. A positive ITPC result would support fixational phase
%   alignment, not the full active sensing account.
%
% Output (when implemented)
%   paths.features/GCP_ms_locked_ITPC.mat
%   figures under paths.figures/exploratory/ms_locked/

%% Setup
startup
[subjects, paths] = setup('GCP');
if ispc
    addpath('C:\Users\Administrator\Documents\GitHub\GCP\_BACKUP\_PILOT');
else
    addpath('/Users/Arne/Documents/GitHub/GCP/_BACKUP/_PILOT');
end
[paths, subjects] = gcp_pilot_paths(paths);
subjects = gcp_subject_inclusion(subjects, paths);

return
