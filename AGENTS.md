# Agent Guidelines for EEGLAB

EEGLAB is a MATLAB/Octave toolbox for EEG/MEG/time-series analysis with both GUI
and command-line workflows. Favor small MATLAB-native patches that reuse current
EEGLAB structure and helpers; avoid new architecture, Python/Node tooling, or
repo-wide cleanup unless explicitly requested.

## Start Here

- Use `.agents/skills/eeglab-matlab-development/SKILL.md` for MATLAB source
  work, `fix-issue` for GitHub issues, `github-pr-review` for reviews, and
  `pull-request` for PR text.
- Read the touched function and nearby analogs first. Search with `rg` before
  adding helpers.
- Key references: `CONTRIBUTING.md`, `eeglab.m`, `functions/adminfunc/eeg_checkset.m`.
- Default branch is `develop`; `CONTRIBUTING.md` says bug fixes target `master`
  and enhancements target `develop`.
- Submodules: `plugins/dipfit`, `plugins/clean_rawdata`, `plugins/ICLabel`,
  `plugins/firfilt`, `plugins/EEG-BIDS`, `tutorial_scripts`. Clone/update with
  `--recurse-submodules` and `git submodule update --init --recursive --remote`.
  Do not edit submodule contents or move submodule pointers casually.

## Run And Validate

Startup modes: `eeglab` GUI, `eeglab nogui` headless, `eeglab redraw`,
`eeglab rebuild`, `eeglab versions`.

```bash
matlab -batch "cd('/path/to/eeglab'); eeglab('nogui'); <commands>"
octave --quiet --eval "cd('/path/to/eeglab'); eeglab('nogui'); <commands>"
```

Smoke check:

```bash
matlab -batch "cd('/path/to/eeglab'); eeglab('nogui'); EEG = pop_loadset('filename','eeglab_data.set','filepath','sample_data/'); EEG = eeg_checkset(EEG);"
```

If MATLAB/Octave, display, license, or data are unavailable, say exactly what was
not run. Do not claim validation from code reading alone.

## Repo Map

- `eeglab.m`: startup, menus, global GUI state, plugin discovery/loading.
- `functions/popfunc/`: `pop_*` GUI/script wrappers and many `eeg_*` functions.
- `functions/adminfunc/`: validation, options, history, dataset store/retrieve,
  plugin admin (`eeg_checkset`, `eeg_eval`, `eeg_store`, `vararg2str`).
- `functions/guifunc/`: `inputgui`, `supergui`, dialogs, channel selection.
- `functions/sigprocfunc/`, `timefreqfunc/`, `statistics/`: core processing.
- `functions/miscfunc/`: import/export, channel/event/ICA/numerical utilities.
- `functions/studyfunc/`: `STUDY`/`ALLEEG` multi-subject analysis.
- `functions/@eegobj`, `@memmapdata`, `@mmo`: MATLAB class folders.
- `plugins/*/eegplugin_*.m`: plugin menu registration.
- `sample_data/`, `sample_data/test_data/`, `sample_locs/`: validation data.

Function categories:
- `pop_*`: user/menu wrappers; no args usually opens GUI; key/value args must
  stay scriptable; return `[EEG, com]` or `[EEG, LASTCOM]` for history.
- `eeg_*`: EEG structure/admin functions (`eeg_checkset`, `eeg_epoch`,
  `eeg_store`, etc.).
- Processing functions: no dialogs; direct algorithms (`runica`, `topoplot`,
  `spectopo`, `eegfilt`, time-frequency functions).

## EEG Structure Invariants

Data is channel-major: continuous `[nbchan x pnts]`, epoched
`[nbchan x pnts x trials]`; `EEG.data` can also reference disk/memmap data.
Event latencies are 1-based sample points, not seconds. Call the relevant
`eeg_checkset` mode after changing data, dimensions, events, epochs, channels,
or ICA fields; use `eeg_checkset(EEG, 'eventconsistency')` after event edits.

| Field group | Key fields |
| --- | --- |
| Dimensions/time | `data`, `nbchan`, `pnts`, `trials`, `srate`, `xmin`, `xmax`, `times` |
| Channels/ref | `chanlocs`, `urchanlocs`, `chaninfo`, `ref`, `splinefile` |
| Events/epochs | `event`, `urevent`, `epoch`, `eventdescription`, `epochdescription` |
| ICA/components | `icaweights`, `icasphere`, `icawinv`, `icaact`, `icachansind`, `dipfit` |
| Rejection/stats | `reject`, `stats`, `specdata`, `specica` |
| STUDY metadata | `subject`, `group`, `condition`, `run`, `session` |
| Misc/save | `etc`, `comments`, `history`, `saved`, `filename`, `filepath` |

When selecting/removing channels, epochs, time ranges, or components, update all
coupled fields. Watch first/last samples, boundary events, empty `icaact`, ICA
rank, `icachansind`, `urevent` links, `EEG.saved`, and on-disk data references.

## MATLAB Patterns To Preserve

- Prefer existing helpers: `finputcheck`, `vararg2str`, `eeg_checkset`,
  `eeg_eval`, `eeg_store`, `eeg_retrieve`, `eeg_decodechan`, `eeg_mergelocs`,
  `fastif`, `inputgui`, `supergui`, `questdlg2`, `pophelp`.
- `pop_*` changes must preserve GUI cancel behavior (`com = ''`), scripted
  behavior, history strings, and multi-dataset paths via `eeg_eval` when present.
- Use `inputgui`/`supergui` for EEGLAB-style dialogs unless the touched file has
  a clear different precedent.
- Keep old compatibility code unless removal is the task. EEGLAB supports old
  datasets, old MATLAB releases, Octave command-line use, and plugins.
- MATLAB style: 2-space indents, spaces not tabs, no space before `(` in calls,
  one space after commas, preserve help/license/history blocks.
- Do not add defensive clutter for states already guaranteed by `finputcheck`,
  callers, or `eeg_checkset`.

## Menu To Function Map

| Menu | Functions |
| --- | --- |
| File | `pop_loadset`, `pop_saveset`, `pop_fileio`, `pop_biosig`, `pop_importdata`, `pop_importevent`, `pop_importepoch`, `pop_export`, `pop_editoptions`, `pop_importbids` |
| Edit | `pop_editset`, `pop_chanedit`, `pop_editeventfield`, `pop_editeventvals`, `pop_select`, `pop_rmdat`, `pop_selectevent`, `pop_mergeset` |
| Tools/preprocess | `pop_resample`, `pop_eegfilt`, `pop_eegfiltnew`, `pop_reref`, `pop_interp`, `pop_eegplot`, `pop_rejchan`, `pop_rejcont`, `pop_autorej`, `pop_runica`, `pop_subcomp`, `pop_epoch`, `pop_rmbase`, `pop_clean_rawdata`, `pop_iclabel`, `pop_icflag` |
| Epoch rejection | `pop_eegthresh`, `pop_rejtrend`, `pop_jointprob`, `pop_rejkurt`, `pop_rejspec` |
| Plot | `pop_eegplot`, `pop_spectopo`, `pop_prop`, `pop_erpimage`, `pop_timtopo`, `pop_topoplot`, `pop_headplot`, `pop_envtopo`, `pop_newtimef` |
| STUDY | `pop_study`, `pop_studywizard`, `pop_studyerp`, `pop_loadstudy`, `pop_savestudy`, `pop_studydesign`, `pop_precomp`, `pop_preclust`, `pop_clust` |

## Plugin Notes

Plugins live under `plugins/` and register through `eegplugin_<name>.m`.
Programmatic install form:

```matlab
plugin_askinstall('ICLabel', 'iclabel', 0);
plugin_askinstall('clean_rawdata', 'clean_artifacts', 0);
plugin_askinstall('firfilt', 'pop_eegfiltnew', 0);
plugin_askinstall('picard', 'picard', 0);
plugin_askinstall('dipfit', 'pop_dipfit_settings', 0);
```

`clean_rawdata` order: flatlines -> high-pass -> bad channels -> ASR bursts ->
bad windows. Any criterion can be `'off'`. Results: `EEG.etc.clean_channel_mask`
and `EEG.etc.clean_sample_mask`.

| clean_rawdata key | Default / note |
| --- | --- |
| `FlatlineCriterion` | `5` seconds |
| `ChannelCriterion` | `0.8`, needs channel locations |
| `LineNoiseCriterion` | `4` SD |
| `Highpass` | `[0.25 0.75]`, use `'off'` if already filtered |
| `BurstCriterion` | `20` GUI conservative; `5` aggressive; `40` mild first pass |
| `BurstRejection` | `'on'` rejects periods; `'off'` corrects via ASR |
| `WindowCriterion` | `0.25` contaminated-channel fraction |
| `WindowCriterionTolerances` | `[-Inf 7]` |
| `Distance` | `'Euclidian'` or `'Riemannian'` |
| `channels_ignore` | labels such as `{'ECG'}` |

`ICLabel`: `pop_iclabel(EEG, 'default')` stores
`EEG.etc.ic_classification.ICLabel.classifications` as `[nComponents x 7]`
probabilities `[Brain Muscle Eye Heart LineNoise ChannelNoise Other]`.
Versions: `'default'`, `'lite'`, `'beta'`. `pop_icflag` uses a `[7 x 2]`
threshold matrix; `NaN NaN` skips a class, for example:

```matlab
EEG = pop_icflag(EEG, [NaN NaN; 0.9 1; 0.9 1; NaN NaN; NaN NaN; NaN NaN; NaN NaN]);
EEG = pop_subcomp(EEG, find(EEG.reject.gcompreject), 0);
```

## Workflow References

Typical ERP pipeline: load/import data -> channel locations -> remove non-EEG
channels -> average reference -> `pop_clean_rawdata` -> re-reference -> ICA
with rank handling -> ICLabel/ICFlag -> remove components -> epoch -> baseline
-> save.

```matlab
EEG = pop_loadset('filename', 'data.set', 'filepath', '/path/');
EEG = pop_chanedit(EEG, 'lookup', 'standard-10-5-cap385.elp');
EEG = pop_select(EEG, 'nochannel', {'EXG1','EXG2','EXG3','ECG','EMG'});
EEG = pop_reref(EEG, []);
EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion', 5, 'ChannelCriterion', 0.8, ...
  'LineNoiseCriterion', 4, 'Highpass', [0.25 0.75], 'BurstCriterion', 20, ...
  'WindowCriterion', 0.25, 'BurstRejection', 'on', 'Distance', 'Euclidian', ...
  'WindowCriterionTolerances', [-Inf 7]);
EEG = pop_reref(EEG, []);
EEG = pop_runica(EEG, 'icatype', 'runica', 'options', {'pca', -1});
EEG = pop_iclabel(EEG, 'default');
EEG = pop_icflag(EEG, [NaN NaN; 0.9 1; 0.9 1; NaN NaN; NaN NaN; NaN NaN; NaN NaN]);
EEG = pop_subcomp(EEG, find(EEG.reject.gcompreject), 0);
EEG = pop_epoch(EEG, {'xxx','yyy'}, [-1 2], 'epochinfo', 'yes');
EEG = pop_rmbase(EEG, [-1000 0]);
EEG = pop_saveset(EEG, 'filename', 'processed.set', 'filepath', '/path/');
```

ICA: `pop_runica` supports `'runica'` (default Infomax), `'picard'` (plugin,
same objective, faster), `'binica'`, `'jader'`, `'sobi'`. Best practice:
high-pass continuous data at 1-2 Hz before ICA, do not baseline-correct before
ICA, use rank reduction after average reference (`'pca', -1` or
`EEG.nbchan - 1`), and train on continuous data when possible.

Filtering: prefer firfilt `pop_eegfiltnew(EEG, 'locutoff', 1)` /
`'hicutoff', 40`; legacy `pop_eegfilt(EEG, 1, 0)` or `(EEG, 0, 40)` still
exists. Filter continuous data before epoching.

Reference/interpolation: `pop_reref(EEG, [])` average reference,
`pop_reref(EEG, [1 2])` indices, `pop_reref(EEG, 'Cz')` label. Average reference
reduces rank by 1. Interpolate removed channels after ICA/component removal:
`pop_interp(EEG, EEG.urchanlocs, 'spherical')`, channel indices, or another
dataset's `chanlocs`.

Events: add/edit events in samples, then `eeg_checkset(EEG, 'eventconsistency')`.
For 100 ms before an event: `EEG.event(end).latency = oldLatency - 0.1*EEG.srate`.
Import events with `pop_importevent(EEG, 'event', file, 'fields', {'latency','type'})`.

STUDY/BIDS outline:

```matlab
[STUDY, ALLEEG] = pop_importbids(filepath, 'eventtype', 'trial_type', ...
  'bidsevent', 'on', 'bidschanloc', 'on', 'studyName', 'MyStudy');
STUDY = std_makedesign(STUDY, ALLEEG, 1, 'name', 'Design1', ...
  'variable1', 'type', 'values1', {'target','standard'}, ...
  'vartype1', 'categorical', 'subjselect', STUDY.subject);
[STUDY, ALLEEG] = std_precomp(STUDY, ALLEEG, {}, 'savetrials', 'on', ...
  'rmicacomps', 'on', 'interp', 'on', 'recompute', 'on', 'erp', 'on');
STUDY = pop_erpparams(STUDY, 'topotime', 350);
STUDY = std_erpplot(STUDY, ALLEEG, 'channels', {ALLEEG(1).chanlocs.labels}, 'design', 1);
```

## Testing Expectations

No central suite exists. Prefer the narrowest reproducible MATLAB/Octave command
using `sample_data/` or `sample_data/test_data/`; plugin-local tests include
`plugins/ICLabel/run_tests.m`.

Validation checklist by change type:
- `pop_*`: GUI path if display is available, scripted path always, returned
  `com`, cancel path, multi-dataset path if present.
- Events/channels/ICA: boundary samples, first/last event, `urevent`, `epoch`,
  channel labels, `nbchan`, `chanlocs`, `icachansind`, empty `icaact`, matrix
  sizes, and `eeg_checkset`.
- GUI: labels/order/tags/callbacks/help match the menu workflow; no GUI-only
  behavior should break command-line use.

## GitHub And Commits

- Never credit yourself or AI tools in commits, comments, or PR text.
- Stage only intended files; this repo may contain unrelated untracked local
  checkouts or tools.
- Use concise commits and terse GitHub comments that state exactly what was
  tested.
- Add `agent-generated` only for repository automation-created issues/PRs, not
  human-directed interactive work.
