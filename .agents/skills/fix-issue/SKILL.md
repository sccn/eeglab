---
name: fix-issue
description: End-to-end workflow for fixing a GitHub issue in sccn/eeglab. Use when asked to fix, investigate, or resolve a GitHub issue.
---

# Skill: Fix GitHub Issue

You are to fix the GitHub issue indicated by the user.

Read first:

@AGENTS.md
@.agents/skills/eeglab-matlab-development/SKILL.md

## Writing Style For GitHub Comments

Be terse. Every sentence must add information.

- No preamble, filler, or generic investigation language.
- Do not repeat the issue body.
- Link or name relevant files instead of narrating obvious code.
- Keep progress comments short and update one comment in place when running from
  automation.

## Research

Use `gh` to fetch the issue. Read enough source to identify the relevant
function, helper chain, plugin, or GUI path.

Before coding, run a duplicate-work preflight:

```bash
gh pr list --state open --search "<issue-id or keyword>"
gh issue list --state open --search "<keyword>"
```

If overlapping work exists, do not open a parallel implementation. Comment with
what you found and either hand off to the existing work or scope your fix to a
non-overlapping follow-up.

Post one comment with these sections before implementing:

```md
# Research

<2-3 sentences explaining the root cause or missing behavior.>

**Relevant code**
- `path/to/file.m#Lx-Ly` - short annotation

# Proposed Fix

<One short paragraph explaining the smallest fix and focused validation.>
```

For larger feature or design work, do not guess. Write a concise design note in
the issue instead of opening a PR immediately. Include the problem, approach,
affected files, tests, tradeoffs, and open questions.

## Implementation

- Create a branch named `agent/{YYYYMMDD}-fix-{issue-id}`.
- Reproduce the issue with a focused MATLAB or Octave command when possible.
- Add a regression test or reproducibility command when practical. EEGLAB has no
  single central test suite, so prefer a minimal scriptable command using
  `sample_data/` or plugin-local tests.
- Make the smallest local MATLAB change.
- Preserve `pop_*` GUI, command-line, and history behavior.
- Keep EEG structure fields consistent and validate with `eeg_checkset`.

## Testing

Use the narrowest relevant validation first:

```bash
matlab -batch "cd('/path/to/eeglab'); eeglab('nogui'); <focused commands>"
```

```bash
octave --quiet --eval "cd('/path/to/eeglab'); eeglab('nogui'); <focused commands>"
```

Run plugin tests when relevant, for example:

```bash
matlab -batch "cd('/path/to/eeglab/plugins/ICLabel'); run_tests"
```

Do not weaken assertions or change expected numerical behavior without proving
the old expectation was wrong.

## Uploading

When validation is complete, open a PR following
`.agents/skills/pull-request/SKILL.md` exactly. Attach it to the issue with a
short comment summarizing the fix and the validation command.

## Verify CI Status

After opening the PR, monitor checks:

```bash
gh pr view <number> --json statusCheckRollup
```

Fix actionable failures before considering the task complete.

## Tasks

- [ ] Fetch issue information.
- [ ] Research relevant EEGLAB source and helper chain.
- [ ] Run duplicate-work preflight.
- [ ] Post Research and Proposed Fix.
- [ ] Create issue branch.
- [ ] Reproduce or define a focused validation command.
- [ ] Implement the smallest fix.
- [ ] Run MATLAB or Octave validation, plus plugin tests if relevant.
- [ ] Push branch.
- [ ] Open PR using the pull-request skill.
- [ ] Monitor CI and address failures.

## Rules

Never credit yourself or AI tools in commits, comments, or PR descriptions.
