# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is EEGLAB?

EEGLAB is an open-source MATLAB toolbox for processing electrophysiological signals (EEG, MEG, and other time series data). It provides a GUI and command-line functions for continuous and event-related data analysis, including independent component analysis (ICA).

## Running MATLAB

```bash
# Non-interactive batch mode (preferred for automation)
/Applications/MATLAB_R2025a.app/bin/matlab -batch "command_here"

# Example: run EEGLAB without GUI and process data
/Applications/MATLAB_R2025a.app/bin/matlab -batch "cd('/Users/arno/GitHub/core_eeg/eeglab'); eeglab nogui; your_script"
```

Startup options: `eeglab` (full GUI), `eeglab nogui` (headless), `eeglab redraw` (refresh GUI), `eeglab rebuild` (close and rebuild).

## Git Workflow

- **`develop`** - Main and default branch
- Submodules: `dipfit`, `clean_rawdata`, `ICLabel`, `firfilt`, `EEG-BIDS`, `tutorial_scripts`
- Clone with `--recurse-submodules`; update with `git submodule update --init --recursive --remote`

## Code Architecture

### The EEG Structure

All processing revolves around the `EEG` struct:

| Field | Type | Description |
|-------|------|-------------|
| `data` | `[chan x pts]` or `[chan x pts x epochs]` | Raw data matrix |
| `nbchan`, `pnts`, `trials` | int | Dimensions (trials=1 for continuous) |
| `srate` | float | Sampling rate in Hz |
| `xmin`, `xmax` | float | Epoch time bounds in seconds |
| `times` | vector | Latency vector in milliseconds |
| `chanlocs` | struct array | Channel names/locations |
| `event` | struct array | Events with `.type`, `.latency`, `.duration` |
| `urevent` | struct array | Original events before any rejection |
| `epoch` | struct array | Epoch metadata (only when epoched) |
| `ref` | string/int | Reference type (`'common'`, `'averef'`, channel index) |
| `icaweights` | matrix | ICA unmixing weights |
| `icasphere` | matrix | ICA sphering matrix |
| `icawinv` | matrix | ICA inverse (mixing) matrix |
| `icaact` | matrix | Component activations (may be empty; recomputed on demand) |
| `dipfit` | struct | Dipole model for ICA components |
| `reject` | struct | Rejection marks (`.gcompreject` = flagged components) |
| `etc` | struct | Miscellaneous (ICLabel results stored here) |
| `history` | cell | Command history for reproducibility |

Always call `eeg_checkset(EEG)` after modifying the structure to validate and recompute derived fields.

### Three Function Categories

1. **`pop_*` functions** (`functions/popfunc/`): GUI wrappers that show dialogs, call processing functions, return `[EEG, LASTCOM]`. Entry points from menus.
2. **`eeg_*` functions** (`functions/adminfunc/`, `functions/popfunc/`): Structure manipulation and validation (`eeg_checkset`, `eeg_epoch`, `eeg_store`).
3. **Processing functions** (`functions/sigprocfunc/`, `functions/timefreqfunc/`): Direct signal processing (`runica`, `topoplot`, `spectopo`, `eegfilt`).

### Plugin Architecture

Plugins live in `plugins/` and register via `eegplugin_[name].m`.

Browse all available plugins: https://sccn.ucsd.edu/eeglab/plugin_uploader/plugin_list_all.php

Install plugins programmatically:
```matlab
% plugin_askinstall(plugin_name, plugin_function, interactive)
% interactive: 0 = install silently, 1 = prompt user
plugin_askinstall('ICLabel', 'iclabel', 0);        % install ICLabel
plugin_askinstall('clean_rawdata', 'clean_artifacts', 0);
plugin_askinstall('firfilt', 'pop_eegfiltnew', 0);
plugin_askinstall('picard', 'picard', 0);
plugin_askinstall('dipfit', 'pop_dipfit_settings', 0);
```

## EEGLAB Menu-to-Function Reference

### File Menu
| Menu Item | Function |
|-----------|----------|
| Load existing dataset | `pop_loadset` |
| Save current dataset(s) | `pop_saveset` |
| Import data from file | `pop_fileio`, `pop_biosig`, `pop_importdata` |
| Import events | `pop_importevent` |
| Import epoch info | `pop_importepoch` |
| Export data to text | `pop_export` |
| Preferences | `pop_editoptions` |
| Import BIDS dataset | `pop_importbids` (EEG-BIDS plugin) |

### Edit Menu
| Menu Item | Function |
|-----------|----------|
| Dataset info | `pop_editset` |
| Channel locations | `pop_chanedit` |
| Event fields | `pop_editeventfield` |
| Event values | `pop_editeventvals` |
| Select data (channels/time) | `pop_select` |
| Select data using events | `pop_rmdat` |
| Select epochs or events | `pop_selectevent` |
| Append datasets | `pop_mergeset` |

### Tools Menu (Preprocessing)
| Menu Item | Function |
|-----------|----------|
| Change sampling rate | `pop_resample` |
| Basic FIR filter (legacy) | `pop_eegfilt` |
| FIR filter (firfilt plugin) | `pop_eegfiltnew` |
| Re-reference | `pop_reref` |
| Interpolate electrodes | `pop_interp` |
| Inspect/reject by eye | `pop_eegplot` |
| Automatic channel rejection | `pop_rejchan` |
| Automatic continuous rejection | `pop_rejcont` |
| Automatic epoch rejection | `pop_autorej` |
| Decompose data by ICA | `pop_runica` |
| Remove components from data | `pop_subcomp` |
| Extract epochs | `pop_epoch` |
| Remove epoch baseline | `pop_rmbase` |
| Clean Rawdata and ASR | `pop_clean_rawdata` (plugin) |
| Classify components ICLabel | `pop_iclabel` (plugin) |
| Flag components as artifacts | `pop_icflag` (plugin) |

### Epoch Rejection Tools
| Menu Item | Function |
|-----------|----------|
| Reject extreme values | `pop_eegthresh` |
| Reject by linear trend/variance | `pop_rejtrend` |
| Reject by probability | `pop_jointprob` |
| Reject by kurtosis | `pop_rejkurt` |
| Reject by spectra | `pop_rejspec` |

### Plot Menu
| Menu Item | Function |
|-----------|----------|
| Channel data (scroll) | `pop_eegplot` |
| Channel spectra and maps | `pop_spectopo` |
| Channel properties | `pop_prop` |
| Channel ERP image | `pop_erpimage` |
| Channel ERPs with scalp maps | `pop_timtopo` |
| ERP map series (2-D) | `pop_topoplot` |
| ERP map series (3-D) | `pop_headplot` |
| Component activations (scroll) | `pop_eegplot` (with ICA data) |
| Component spectra and maps | `pop_spectopo` (with ICA) |
| Component maps (2-D) | `pop_topoplot` (with ICA) |
| Component properties | `pop_prop` (with ICA) |
| Component ERPs | `pop_envtopo` |
| Time-frequency | `pop_newtimef` |

### STUDY Menu (Multi-Subject)
| Menu Item | Function |
|-----------|----------|
| Create STUDY (loaded datasets) | `pop_study` |
| Browse for datasets | `pop_studywizard` |
| Simple ERP STUDY | `pop_studyerp` |
| Load/Save STUDY | `pop_loadstudy` / `pop_savestudy` |
| Edit STUDY design | `pop_studydesign` |
| Pre-compute statistics | `pop_precomp` |
| Pre-cluster components | `pop_preclust` |
| Cluster components | `pop_clust` |

## Standard Preprocessing Pipeline

Recommended order for ERP analysis:

```matlab
% 1. Load data
EEG = pop_loadset('filename', 'data.set', 'filepath', '/path/');
% or: EEG = pop_fileio('/path/to/data.edf');
% or: [STUDY, ALLEEG] = pop_importbids(bidspath, 'studyName', 'MyStudy');

% 2. Import channel locations (if not already present)
EEG = pop_chanedit(EEG, 'lookup', 'standard-10-5-cap385.elp');

% 3. Remove non-EEG channels (EMG, EOG, ECG, GSR, etc.)
EEG = pop_select(EEG, 'nochannel', {'EXG1','EXG2','EXG3','ECG','EMG'});

% 4. Average reference (before artifact cleaning)
EEG = pop_reref(EEG, []);

% 5. Clean data: remove bad channels, reject bad segments (clean_rawdata)
EEG = pop_clean_rawdata(EEG, ...
    'FlatlineCriterion', 5, ...
    'ChannelCriterion', 0.8, ...
    'LineNoiseCriterion', 4, ...
    'Highpass', [0.25 0.75], ...
    'BurstCriterion', 20, ...
    'WindowCriterion', 0.25, ...
    'BurstRejection', 'on', ...
    'Distance', 'Euclidian', ...
    'WindowCriterionTolerances', [-Inf 7]);

% 6. Re-reference again (after bad channel removal)
EEG = pop_reref(EEG, []);

% 7. Run ICA (pca -1 = auto-reduce for rank-deficient data)
% 'pca', -1 indicate to reduce the dimension by 1 to account for rank decrease by average reference in 6
EEG = pop_runica(EEG, 'icatype', 'runica', 'options', {'pca', -1});

% 8. Classify and flag artifact components (ICLabel)
EEG = pop_iclabel(EEG, 'default');
EEG = pop_icflag(EEG, [NaN NaN; 0.9 1; 0.9 1; NaN NaN; NaN NaN; NaN NaN; NaN NaN]);
%                       Brain   Muscle  Eye    Heart  LineNoise ChanNoise Other

% 9. Remove flagged components
EEG = pop_subcomp(EEG, find(EEG.reject.gcompreject), 0);

% 10. Extract epochs
EEG = pop_epoch(EEG, {'xxx','yyy'}, [-1 2], 'epochinfo', 'yes');
EEG = eeg_checkset(EEG);

% 11. Remove baseline
EEG = pop_rmbase(EEG, [-1000 0]);

% 12. Save
EEG = pop_saveset(EEG, 'filename', 'processed.set', 'filepath', '/path/');
```

## Key Plugin Reference: clean_rawdata

Automated artifact rejection on continuous data via Artifact Subspace Reconstruction (ASR).

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `FlatlineCriterion` | 5 | Max flatline duration (seconds) before channel removal |
| `ChannelCriterion` | 0.8 | Min correlation with neighbors (0-1); requires channel locations |
| `LineNoiseCriterion` | 4 | Max line noise relative to population (std devs) |
| `Highpass` | [0.25 0.75] | Transition band for ~0.5 Hz high-pass. Use `'off'` if already filtered |
| `BurstCriterion` | 20 (GUI) | ASR threshold (std devs). 5=aggressive, 20=conservative |
| `BurstRejection` | `'on'` | `'on'`=reject segments, `'off'`=correct via ASR (preserves data length) |
| `WindowCriterion` | 0.25 | Max fraction of contaminated channels per window |
| `WindowCriterionTolerances` | [-Inf 7] | Power tolerance bounds for window criterion |
| `Distance` | `'Euclidian'` | `'Euclidian'` or `'Riemannian'` distance metric |
| `channels_ignore` | [] | Cell array of channel labels to exclude (e.g., `{'ECG'}`) |

Any parameter can be set to `'off'` to skip that step. Internal processing order: flatlines -> highpass -> bad channels -> ASR bursts -> bad windows.

Results stored in `EEG.etc.clean_channel_mask` and `EEG.etc.clean_sample_mask`.

### Two-Pass Strategy if data very noisy

For best ICA quality, use two passes:
1. **Pass 1 (mild):** `'BurstCriterion', 40` - remove only bad channels and extreme artifacts
2. **Run ICA + ICLabel** on mildly cleaned data
3. **Pass 2 (aggressive):** `'BurstCriterion', 20` - clean remaining artifacts

## Key Plugin Reference: ICLabel

Deep-learning classifier for ICA component labeling. Trained on >500,000 crowd-sourced labeled components.

### Classification

```matlab
EEG = pop_iclabel(EEG, 'default');
% Results in: EEG.etc.ic_classification.ICLabel.classifications  (N_components x 7 matrix)
% Columns:    [Brain, Muscle, Eye, Heart, LineNoise, ChannelNoise, Other]
% Each row sums to 1.0
```

Versions: `'default'` (recommended), `'lite'` (faster, no autocorrelation), `'beta'` (legacy).

### Flagging and Removal

`pop_icflag` threshold matrix is `[7x2]`: each row = `[min max]` probability for `[Brain, Muscle, Eye, Heart, LineNoise, ChannelNoise, Other]`. Use `NaN NaN` to skip a category:

```matlab
% Flag Muscle (>90%) and Eye (>90%) artifacts
EEG = pop_icflag(EEG, [NaN NaN; 0.9 1; 0.9 1; NaN NaN; NaN NaN; NaN NaN; NaN NaN]);

% Flag anything with <20% Brain probability
EEG = pop_icflag(EEG, [0 0.2; NaN NaN; NaN NaN; NaN NaN; NaN NaN; NaN NaN; NaN NaN]);

% Remove flagged components
EEG = pop_subcomp(EEG, find(EEG.reject.gcompreject), 0);
```

## ICA Reference

### Algorithms Available via pop_runica

| Algorithm | `'icatype'` value | Notes |
|-----------|-------------------|-------|
| Infomax | `'runica'` | Default MATLAB implementation |
| Picard | `'picard'` | Faster convergence, same objective as runica. Requires separate plugin. Recommended. |
| Binary Infomax | `'binica'` | Compiled C, faster than runica |
| JADE | `'jader'` | |
| SOBI | `'sobi'` | Second-order blind identification |

### ICA Best Practices

- **High-pass filter at 1-2 Hz** before ICA (critical for decomposition quality)
- Clean_rawdata's 0.5 Hz default is a compromise; use `pop_eegfiltnew(EEG, 'locutoff', 1)` for tighter
- Average referencing reduces rank by 1: use `'pca', -1` (auto) or `'pca', EEG.nbchan - 1`
- Run on continuous data (not epoched) for maximum training samples
- Do NOT baseline-correct before ICA
- Data requirement: roughly >30*N^2 samples where N = number of channels

## Re-Referencing

```matlab
EEG = pop_reref(EEG, []);              % Average reference
EEG = pop_reref(EEG, [1 2]);           % Reference to channels 1 and 2
EEG = pop_reref(EEG, 'Cz');            % Reference to named channel
```

Average reference reduces data rank by 1 (important for ICA dimensionality). Re-reference BEFORE ICA. Multiple average references cancel earlier ones.

## Channel Interpolation

Interpolate removed or bad channels using spherical spline (default) or other methods. Best done **after ICA component removal** -- interpolating before ICA introduces artificial data that degrades decomposition quality.

```matlab
% Interpolate missing channels from a full-montage reference
% (urchanlocs preserves the original channel list before any were removed)
EEG = pop_interp(EEG, EEG.urchanlocs, 'spherical');

% Interpolate specific channels by index
EEG = pop_interp(EEG, [12 48], 'spherical');

% Interpolate using a different dataset's channel locations as template
EEG = pop_interp(EEG, ALLEEG(1).chanlocs, 'spherical');
```

Typical pipeline position: clean_rawdata (removes bad channels) -> re-reference -> ICA -> ICLabel -> remove components -> **interpolate** -> re-reference (again, optional) -> epoch.

## Filtering

```matlab
% FIR filter (firfilt plugin - preferred)
EEG = pop_eegfiltnew(EEG, 'locutoff', 1);           % 1 Hz high-pass
EEG = pop_eegfiltnew(EEG, 'hicutoff', 40);          % 40 Hz low-pass
EEG = pop_eegfiltnew(EEG, 'locutoff', 1, 'hicutoff', 40);  % Bandpass

% Legacy FIR filter
EEG = pop_eegfilt(EEG, 1, 0);    % 1 Hz high-pass
EEG = pop_eegfilt(EEG, 0, 40);   % 40 Hz low-pass
```

Always filter continuous data before epoching.

## Event Manipulation

```matlab
% Add new events programmatically
for i = 1:length(EEG.event)
    if strcmpi(EEG.event(i).type, 'stimulus')
        EEG.event(end+1) = EEG.event(i);
        EEG.event(end).latency = EEG.event(i).latency - 0.1*EEG.srate;  % 100ms before
        EEG.event(end).type = 'cue';
    end
end
EEG = eeg_checkset(EEG, 'eventconsistency');

% Import events from file
EEG = pop_importevent(EEG, 'event', 'events.txt', 'fields', {'latency','type'});
```

Event latencies are in sample points (1-indexed). Convert to seconds: `latency_sec = EEG.event(i).latency / EEG.srate`.

## STUDY-Level Analysis (Multi-Subject)

```matlab
% Create STUDY from BIDS
[STUDY, ALLEEG] = pop_importbids(filepath, 'eventtype', 'trial_type', ...
    'bidsevent', 'on', 'bidschanloc', 'on', 'studyName', 'MyStudy');

% Create STUDY design
STUDY = std_makedesign(STUDY, ALLEEG, 1, 'name', 'Design1', ...
    'variable1', 'type', 'values1', {'target','standard'}, ...
    'vartype1', 'categorical', 'subjselect', STUDY.subject);

% Precompute measures
[STUDY, ALLEEG] = std_precomp(STUDY, ALLEEG, {}, 'savetrials', 'on', ...
    'rmicacomps', 'on', 'interp', 'on', 'recompute', 'on', 'erp', 'on');

% Plot
STUDY = pop_erpparams(STUDY, 'topotime', 350);
STUDY = std_erpplot(STUDY, ALLEEG, 'channels', {ALLEEG(1).chanlocs.labels}, 'design', 1);
```

## MATLAB Code Style (EEGLAB conventions)

- 2-space indentation (spaces, not tabs)
- No space between function name and parenthesis: `eeg_checkset(EEG)` not `eeg_checkset (EEG)`
- One space after commas in argument lists
- `pop_*` functions return `[EEG, LASTCOM]` where LASTCOM is the command string for history

## Testing

No centralized test suite. Testing is per-plugin (`plugins/ICLabel/run_tests.m`, etc.) and manual:
1. `eeglab nogui` loads without errors
2. Load sample data: `EEG = pop_loadset('filename', 'eeglab_data.set', 'filepath', 'sample_data/')`
3. Run the modified function
4. Validate with `EEG = eeg_checkset(EEG)`
