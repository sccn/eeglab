---
name: pull-request
description: Authoring pull requests for sccn/eeglab. Use when creating or updating a PR, and whenever changing a branch that is already associated with a PR.
---

# Skill: Author a Pull Request

This skill defines the PR format and pre-push expectations for EEGLAB work.
Follow it when creating or updating a PR.

Read first:

@AGENTS.md

## PR Description Format

The PR description should be plain text and short.

Title: Short imperative sentence. Optional scope tag in brackets.

Body: 1-3 sentences stating what changed and why. End with an issue link if one
exists.

Hard rules:

- No markdown headers, tables, images, or checkboxes.
- No `## Summary`, `## Test plan`, or section labels.
- No filler phrases such as "This PR" or "Summary of changes".
- No emoji.
- Under about 80 words total.
- Never credit yourself or AI tools.

Example:

```text
Title: [EEG] Fix event latency update after selection

Body:
Keep event latencies consistent when selecting a trailing continuous-data
window. Adds a sample-data regression command covering first and last event
latencies.

Fixes #1234
```

## Issue Linking

Use `Fixes #NNNN` only when the change resolves a pre-existing issue. Use
`Part of #NNNN` for partial work. Omit issue links for user-directed work with
no issue.

## Pre-Push Checklist

Run the narrowest meaningful validation before pushing:

- Core smoke:

  ```bash
  matlab -batch "cd('/path/to/eeglab'); eeglab('nogui'); EEG = pop_loadset('filename', 'eeglab_data.set', 'filepath', 'sample_data/'); EEG = eeg_checkset(EEG);"
  ```

- Octave smoke when MATLAB is unavailable:

  ```bash
  octave --quiet --eval "cd('/path/to/eeglab'); eeglab('nogui'); EEG = pop_loadset('filename', 'eeglab_data.set', 'filepath', 'sample_data/'); EEG = eeg_checkset(EEG);"
  ```

- Focused test or command for the function, plugin, import path, GUI path, or
  STUDY path that changed.

If validation cannot run because MATLAB, Octave, a display, a license, or test
data are unavailable, say that clearly in the PR body or final status.

## Creating The PR

Unless the user says otherwise, push directly to a branch on `sccn/eeglab` when
permissions allow. Do not push to a fork unless direct push is unavailable or
requested.

Use:

```bash
gh pr create \
  --title "<title>" \
  --body "<plain text body>"
```

Add the `agent-generated` label only when repository automation creates the PR.
Do not add it when a human asks an interactive agent to open or update a PR.

## See Also

- `.agents/skills/fix-issue/`
- `.agents/skills/github-pr-review/`
- `AGENTS.md`
