---
name: eeglab-matlab-development
description: Use when implementing, fixing, or reviewing EEGLAB MATLAB source, including pop_ functions, EEG structure updates, GUI dialogs, plugin registration, STUDY functions, event/channel/ICA logic, or MATLAB/Octave validation.
---

# EEGLAB MATLAB Development

Use this skill for source changes in the EEGLAB MATLAB repository. The goal is
to make the smallest correct MATLAB-native change that fits the existing
toolbox conventions.

## Start From Existing Code

Read first:

@AGENTS.md

Then inspect the closest existing implementation before editing. Good starting
points:

- `eeglab.m` for startup, menus, redraw behavior, plugin loading, and global GUI
  state.
- `functions/popfunc/pop_*.m` for user-facing wrappers and command history.
- `functions/adminfunc/eeg_checkset.m` for EEG structure invariants.
- `functions/adminfunc/eeg_eval.m` for multi-dataset execution.
- `functions/guifunc/inputgui.m` and `functions/guifunc/supergui.m` for dialog
  construction.
- `plugins/*/eegplugin_*.m` for plugin menu registration.
- The neighboring function in the same folder for style, parsing, and output
  conventions.

Before adding a helper, search for an existing one:

```bash
rg "function .*<name>" functions plugins
rg "<option-or-field-name>" functions plugins
```

## Implementation Rules

- Keep changes small and local. Patch the branch that is wrong instead of
  rewriting mature functions.
- Prefer `finputcheck`, `vararg2str`, `eeg_checkset`, `eeg_eval`,
  `eeg_decodechan`, `eeg_store`, `eeg_retrieve`, `inputgui`, and `supergui`
  over new local frameworks.
- Preserve both GUI and command-line usage for `pop_*` functions.
- Preserve returned command strings. Build them from actual options, usually
  with `vararg2str`.
- Do not add Python, Node, generated project files, or non-MATLAB tooling for
  normal EEGLAB code changes.
- Be careful with plugins: `plugins/ICLabel`, `plugins/clean_rawdata`,
  `plugins/dipfit`, `plugins/firfilt`, `plugins/EEG-BIDS`, and
  `tutorial_scripts` are submodules in the superproject.
- Avoid new `eval` unless matching an established EEGLAB pattern is necessary.
- Keep 1-based MATLAB indexing and 1-based event latency semantics explicit.

## Common Change Patterns

### pop_ wrapper change

1. Read the full `pop_<name>.m` function and any called processing function.
2. Confirm the no-argument GUI path and the scripted key/value path.
3. Keep `com = ''` behavior on cancel or early return.
4. Use `finputcheck` for new options when nearby code does.
5. If multiple datasets are supported, keep or add the `eeg_eval` path.
6. Validate with a command-line sample-data call and `eeg_checkset`.

### EEG structure change

1. Identify every field coupled to the modified data, such as `nbchan`,
   `pnts`, `trials`, `times`, `event`, `urevent`, `epoch`, `chanlocs`,
   `icachansind`, `icaact`, `icaweights`, `icasphere`, and `icawinv`.
2. Update related fields together.
3. Use `eeg_checkset(EEG)` or the specific mode already used nearby.
4. Test boundary cases: first/last sample, empty events, single channel, epoched
   versus continuous data, and empty `icaact` where relevant.

### GUI change

1. Use existing `inputgui` and `supergui` patterns unless the touched file
   already uses another approach.
2. Preserve labels, callbacks, tags, button names, help commands, and command
   history behavior.
3. Keep the scriptable path independent from the GUI path.
4. Validate the no-argument GUI manually in MATLAB when a display is available.

### Plugin change

1. Check whether the plugin is a submodule and inspect its own status.
2. Read `eegplugin_<name>.m`, the relevant `pop_*` wrapper, and lower-level
   processing functions.
3. Keep plugin menu registration and dependency checks consistent with sibling
   plugin files.
4. Run plugin tests if present.

## Validation Commands

Prefer MATLAB when available for core EEGLAB behavior:

```bash
matlab -batch "cd('/path/to/eeglab'); eeglab('nogui'); EEG = pop_loadset('filename', 'eeglab_data.set', 'filepath', 'sample_data/'); EEG = eeg_checkset(EEG);"
```

Use Octave for command-line smoke coverage when MATLAB is unavailable:

```bash
octave --quiet --eval "cd('/path/to/eeglab'); eeglab('nogui'); EEG = pop_loadset('filename', 'eeglab_data.set', 'filepath', 'sample_data/'); EEG = eeg_checkset(EEG);"
```

For a changed function, run a focused command that exercises that function on
`sample_data/eeglab_data.set` or a relevant file under `sample_data/test_data/`.

Examples:

```bash
matlab -batch "cd('/path/to/eeglab'); eeglab('nogui'); EEG = pop_loadset('filename', 'eeglab_data.set', 'filepath', 'sample_data/'); EEG = pop_select(EEG, 'channel', 1:4); EEG = eeg_checkset(EEG); assert(EEG.nbchan == 4);"
```

```bash
matlab -batch "cd('/path/to/eeglab/plugins/ICLabel'); run_tests"
```

If a GUI cannot be validated in batch mode, say so and still validate the
non-GUI command path.

## Checklist

- [ ] Read `AGENTS.md`.
- [ ] Found the closest existing implementation.
- [ ] Chose the smallest local patch.
- [ ] Preserved `pop_*` command-line, GUI, and history behavior where relevant.
- [ ] Kept EEG structure fields consistent.
- [ ] Considered MATLAB and Octave compatibility.
- [ ] Ran focused validation or stated why it could not run.
