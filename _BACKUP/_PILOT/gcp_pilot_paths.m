function [paths, subjects] = gcp_pilot_paths(paths)
% Point a setup('GCP') path struct at the archived pilot tree.
% Data and figures stay under Students/Arne/GCP/_PILOT.
% Raw recordings stay in OCC/GCP/data_pilot.

if ispc
    pilotRoot = 'W:\Students\Arne\GCP\_PILOT';
    paths.raw_occ = 'V:\OCC\GCP\data_pilot';
    paths.code = 'C:\Users\Administrator\Documents\GitHub\GCP\_BACKUP\_PILOT';
else
    pilotRoot = '/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/_PILOT';
    paths.raw_occ = '/Volumes/g_psyplafor_methlab_data$/OCC/GCP/data_pilot';
    paths.code = '/Users/Arne/Documents/GitHub/GCP/_BACKUP/_PILOT';
end

paths.root = pilotRoot;
paths.data = fullfile(pilotRoot, 'data');
paths.features = fullfile(paths.data, 'features');
paths.figures = fullfile(pilotRoot, 'figures');
paths.stats = fullfile(paths.data, 'stats');
paths.controls = fullfile(paths.data, 'controls');
paths.logs = fullfile(paths.controls, 'logs');
paths.automagic = fullfile(paths.data, 'automagic');
paths.merged = fullfile(paths.data, 'merged');

subjects = {};
if isfolder(paths.features)
    listing = dir(paths.features);
    names = {listing([listing.isdir]).name};
    names = names(~ismember(names, {'.', '..'}));
    subjects = sort(names);
end
end
