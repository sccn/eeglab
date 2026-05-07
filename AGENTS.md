# Agent Guidelines for EEGLAB

EEGLAB is a MATLAB-native toolbox for electrophysiological signal processing.
Most users interact with it through MATLAB, Octave command-line mode, the EEGLAB
GUI, and EEGLAB plugin conventions. Agent work in this repository should favor
small, direct MATLAB changes that fit the existing structure over new
abstractions, broad rewrites, or cross-language tooling.

Primary references in this repository:
- `README.md` for installation and submodule expectations.
- `CONTRIBUTING.md` for branch, commit, and MATLAB style guidance.
- `eeglab.m` for startup, menus, plugin loading, and user-facing workflow.
- `functions/adminfunc/eeg_checkset.m` for the canonical EEG structure fields.
- `.agents/skills/` for task-specific agent workflows.

## Repo Map

- `eeglab.m`: main EEGLAB entry point, menu construction, startup modes, plugin
  loading, and GUI redraw behavior.
- `functions/popfunc/`: user-facing `pop_*` wrappers. These usually support a
  GUI path when called with only `EEG`, parse key/value arguments for scripted
  use, call lower-level processing functions, and return a command string.
- `functions/adminfunc/`: EEG structure validation, history, options, dataset
  storage, plugin management, and administrative helpers such as
  `eeg_checkset`, `eeg_store`, `eeg_retrieve`, `eeg_eval`, and `vararg2str`.
- `functions/guifunc/`: EEGLAB GUI helpers such as `inputgui`, `supergui`, and
  channel-selection dialogs. Keep GUI construction in these existing patterns.
- `functions/sigprocfunc/`: core signal-processing functions such as ICA,
  filtering, plotting, spectrum, rejection, and interpolation helpers.
- `functions/timefreqfunc/`: time-frequency analysis and plotting functions.
- `functions/miscfunc/`: import/export, channel, event, ICA, and numerical
  utility functions used across the toolbox.
- `functions/statistics/`: statistical helpers and older test scripts.
- `functions/studyfunc/`: STUDY-level functions operating on `STUDY` and
  `ALLEEG` across subjects, conditions, and designs.
- `functions/@eegobj/`, `functions/@memmapdata/`, `functions/@mmo/`: MATLAB
  class-style folders. Preserve MATLAB method dispatch conventions here.
- `plugins/`: bundled submodules. Each plugin registers through
  `eegplugin_<name>.m`. Treat plugin changes as work inside a submodule unless
  the task explicitly targets the superproject pointer or plugin integration.
- `sample_data/` and `sample_data/test_data/`: checked-in data for smoke tests,
  import tests, examples, and reproducibility checks.
- `sample_locs/`: standard channel location files.
- `tutorial_scripts/`: tutorial script submodule.
- `.github/workflows/`: CI and agent automation.

## Before Coding

- Check whether a matching skill exists in `.agents/skills/`. Use
  `.agents/skills/eeglab-matlab-development/SKILL.md` for MATLAB source changes,
  `fix-issue` for GitHub issues, `github-pr-review` for PR review, and
  `pull-request` when preparing PR text.
- Read the existing function, nearby related functions, and the called helper
  chain before editing. EEGLAB has many established helper APIs; prefer them.
- Search before adding code. Look for existing implementations with `rg`, and
  check whether the pattern already exists in a nearby `pop_*`, `eeg_*`, plugin,
  GUI, or STUDY function.
- State assumptions before implementation when a request can mean more than one
  thing. Ask only when a reasonable local assumption would be risky.
- Define the smallest verifiable result. For a bug, reproduce it with a focused
  MATLAB or Octave command if possible, then fix it. For a feature, identify the
  function, menu path, command-line call, and minimal data path to validate.
- Do not introduce Python, Node, package managers, or generated build systems
  for ordinary EEGLAB development. This is a MATLAB and Octave repository.

## Development Principles

- Make the smallest change that solves the request. Avoid broad formatting
  diffs, speculative refactors, compatibility layers, or rewrites of old code
  unless the user explicitly asks for them.
- Prefer existing EEGLAB helpers over new utilities. Common helpers include
  `finputcheck`, `vararg2str`, `eeg_checkset`, `eeg_store`, `eeg_retrieve`,
  `eeg_eval`, `eeg_decodechan`, `eeg_mergelocs`, `fastif`, `inputgui`,
  `supergui`, `questdlg2`, and `pophelp`.
- Keep command-line behavior and GUI behavior aligned. A `pop_*` function should
  remain scriptable without opening a dialog when arguments are supplied.
- Preserve history behavior. `pop_*` functions normally return `com` or
  `LASTCOM` strings built from the actual options used, often via
  `vararg2str`.
- Prefer direct procedural MATLAB code over new classes or deep helper layers.
  Many EEGLAB files are old and long; fit the local style instead of imposing a
  new architecture.
- Avoid hidden global state. When existing EEGLAB globals or options are used,
  keep the use narrow and consistent with nearby code.
- Keep old compatibility code unless the task is specifically to remove it.
  EEGLAB supports old datasets, old MATLAB releases, Octave command-line use,
  and many plugin workflows.
- Do not edit submodules casually. If work is inside `plugins/ICLabel`,
  `plugins/clean_rawdata`, `plugins/dipfit`, `plugins/firfilt`,
  `plugins/EEG-BIDS`, or `tutorial_scripts`, check submodule status and keep the
  superproject pointer changes intentional.

## EEG Structure

Most processing revolves around the `EEG` struct:

- `data`: channel-major numeric data, usually `[nbchan x pnts]` for continuous
  data or `[nbchan x pnts x trials]` for epoched data. It may also be a filename
  or memory-mapped reference when data are stored on disk.
- `nbchan`, `pnts`, `trials`, `srate`, `xmin`, `xmax`, `times`: core dimension
  and timing fields.
- `chanlocs`, `urchanlocs`, `chaninfo`, `ref`: channel labels, locations,
  original channels, and reference metadata.
- `event`, `urevent`, `epoch`, `eventdescription`, `epochdescription`: event
  and epoch metadata. Event latencies are in 1-based sample points.
- `icaweights`, `icasphere`, `icawinv`, `icaact`, `icachansind`, `dipfit`: ICA
  and component metadata. `icaact` may be empty and computed on demand.
- `reject`, `stats`, `specdata`, `specica`: rejection and statistics fields.
- `etc`, `comments`, `history`, `saved`, `filename`, `filepath`: miscellaneous,
  history, save-state, and file location metadata.

After modifying data, dimensions, events, epochs, channel locations, or ICA
fields, call the relevant `eeg_checkset` mode. Use `eeg_checkset(EEG)` for broad
validation, `eeg_checkset(EEG, 'eventconsistency')` after event edits, and
specialized modes only when the existing code path already uses them.

## Function Patterns

- `pop_*` functions usually start with a help block, initialize `com = ''`, show
  help and return when required input is missing, open a GUI only when arguments
  are absent, parse key/value options with `finputcheck`, process multiple
  datasets with `eeg_eval`, update `EEG`, then return a command string.
- Lower-level processing functions should not open dialogs. Keep GUI-specific
  code in `pop_*` wrappers or `functions/guifunc/`.
- GUI code should use `inputgui` or existing EEGLAB GUI helpers unless there is
  a clear local precedent for raw `uicontrol` or figure construction.
- Plugin registration belongs in `eegplugin_<name>.m`. Menu labels, callbacks,
  and version checks should follow existing plugin files.
- STUDY changes should preserve `STUDY`, `ALLEEG`, current design, and saved
  precompute conventions. Inspect nearby `std_*` and `pop_*study*` functions
  before touching STUDY logic.
- Import/export changes should be tested with `sample_data/test_data/` where
  possible and should avoid breaking BIDS or plugin import paths.

## MATLAB Style

- Follow `CONTRIBUTING.md`: 2-space indentation with spaces, no space between a
  function name and `(`, one space after commas, and vertical alignment when it
  genuinely improves readability.
- Preserve existing help-block style, author/history notes, and license header
  conventions in edited files.
- Use `fprintf`, `disp`, `warning`, and `error` consistently with nearby code.
  Avoid noisy output on normal successful paths unless the function already
  reports progress.
- Keep comments useful and sparse. Add comments for non-obvious EEG semantics,
  MATLAB/Octave compatibility, or historical compatibility constraints.
- Avoid new `eval` usage unless it matches an established EEGLAB GUI/history
  pattern and there is no simple safer alternative.
- Avoid vectorization that makes the code harder to audit. Simple loops are
  acceptable when they match the surrounding code and preserve behavior.

## Testing And Validation

There is no single repository-wide test suite. Use the narrowest relevant
MATLAB or Octave validation and broaden when risk requires.

Useful smoke checks:

```bash
matlab -batch "cd('/path/to/eeglab'); eeglab('nogui'); EEG = pop_loadset('filename', 'eeglab_data.set', 'filepath', 'sample_data/'); EEG = eeg_checkset(EEG);"
```

```bash
octave --quiet --eval "cd('/path/to/eeglab'); eeglab('nogui'); EEG = pop_loadset('filename', 'eeglab_data.set', 'filepath', 'sample_data/'); EEG = eeg_checkset(EEG);"
```

Focused validation examples:

- Changed a `pop_*` wrapper: run the command-line call with `sample_data`, check
  returned `EEG`, `com`, and `eeg_checkset`.
- Changed event logic: validate first and last event latencies, boundary events,
  `urevent` links when relevant, and `eeg_checkset(EEG, 'eventconsistency')`.
- Changed channel logic: validate channel labels, `nbchan`, `chanlocs`,
  `urchanlocs`, ICA channel indices, and rejection masks.
- Changed ICA/component logic: validate ICA matrix dimensions, empty `icaact`,
  component removal, and ICLabel or dipfit metadata when relevant.
- Changed GUI code: verify the no-argument GUI path manually in MATLAB when a
  display is available, and separately verify the non-GUI command path.
- Changed a plugin: run that plugin's local tests if present, such as
  `plugins/ICLabel/run_tests.m`, and validate the plugin registration path.

If MATLAB or Octave is unavailable, state that explicitly and run whatever
static or code-reading checks are still meaningful. Do not claim validation that
was not run.

## GitHub, Commits, And Communication

- Never credit yourself or AI tools in commits, comments, or PR descriptions.
- Keep commits scoped to one logical change with a concise message.
- Do not stage unrelated user changes. This repository may contain untracked
  local checkouts or tool directories.
- Use `gh` with narrow JSON fields or explicit flags when inspecting issues and
  PRs. Avoid noisy generic views when a targeted query is enough.
- Agent comments on GitHub issues or PRs should be terse and should state
  exactly what was tested.
- Add the `agent-generated` label only when repository automation creates the PR
  or issue. Do not add it when a human asks an interactive agent to work.

## Failure Patterns To Avoid

- Rewriting a mature MATLAB function instead of patching the relevant branch.
- Introducing a new helper for one call site.
- Adding defensive checks for states that `finputcheck`, `eeg_checkset`, or an
  existing caller already guarantees.
- Breaking the command-line path while changing a GUI dialog.
- Forgetting 1-based event latency semantics.
- Updating data dimensions without keeping `nbchan`, `pnts`, `trials`, `times`,
  events, epochs, ICA fields, or `saved` consistent.
- Editing a plugin submodule without noticing that the superproject pointer will
  need an intentional update.
