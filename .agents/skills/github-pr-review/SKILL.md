---
name: github-pr-review
description: Review a pull request for correctness, MATLAB/Octave behavior, EEGLAB data-structure consistency, GUI/history regressions, and repository convention compliance in sccn/eeglab.
allowed-tools: Bash(gh issue view:*), Bash(gh search:*), Bash(gh issue list:*), Bash(gh pr comment:*), Bash(gh pr diff:*), Bash(gh pr view:*), Bash(gh pr list:*), Bash(gh pr review:*), Bash(git rev-parse:*), mcp__github_inline_comment__create_inline_comment
---

Provide a code review for the given pull request.

You are reviewing EEGLAB, a MATLAB-native toolbox for EEG and related
electrophysiological signal processing.

Read first:

@AGENTS.md

## Review Goals

Prioritize:

1. Correctness bugs and behavioral regressions in changed MATLAB code.
2. EEG structure consistency, including dimensions, events, epochs, channel
   locations, ICA fields, STUDY fields, and save/history fields.
3. GUI and command-line parity for `pop_*` functions.
4. MATLAB and Octave compatibility where the changed code claims or implies it.
5. Focused tests or reproducibility commands for changed behavior.
6. Simple code that follows existing EEGLAB conventions.

## High-Signal Threshold

Flag only issues that are specific, actionable, and tied to changed code.

Flag:

- Code that will fail to parse or run in a common MATLAB or Octave path.
- Incorrect dimensions, indexing, event latencies, channel metadata, ICA fields,
  or STUDY state after the changed operation.
- A broken `pop_*` command-line path, GUI path, or returned command string.
- Missing validation for a realistic regression introduced by the change.
- Unsafe file, path, download, or plugin handling in changed code.
- Clear violations of `AGENTS.md` or scoped directory instructions.

Do not flag:

- Pre-existing issues outside the PR.
- Generic style preferences.
- Broad refactors unrelated to the diff.
- Hypothetical edge cases without a realistic EEGLAB path.
- Issues a normal linter or MATLAB parser would report unless they change
  behavior or block execution.

## Review Workflow

1. Check whether the PR is closed or draft. If so, stop unless the user
   explicitly requested review anyway.
2. Read the PR title, body, changed files, and diff:

   ```bash
   gh pr view <PR> --json title,body,files,state,isDraft,author
   gh pr diff <PR>
   ```

3. Find relevant instructions:

   ```bash
   find .. -name AGENTS.md -o -name CLAUDE.md
   ```

   Apply only instructions in the changed file's directory or parents.

4. Review changed files against the EEGLAB rules in `AGENTS.md`.
5. Validate candidate findings by reading the surrounding source, not just the
   diff.
6. Produce a concise terminal summary. If `--comment` was not provided, stop.
7. If `--comment` was provided, post or update a single PR comment.

## Review Comment Format

Use this format:

```md
## Code review

- Overall assessment: <safe to merge / needs changes / needs more context>
- Highest-risk area: <area>
- Merge recommendation: <safe to merge / needs changes / needs more context>

## Blocking

<findings or None.>

## Important

<findings or None.>

## Nits

<findings or None.>

## Test gaps

<findings or None.>

## EEGLAB notes

<findings or None.>
```

Rules:

- If a section has no findings, write `None.`
- Each finding must include severity, file/function reference, why it matters,
  and a concrete suggested fix.
- Keep the review concise.
- Do not include generic praise.
- If no issues were found, say:

```md
## Code review

- Overall assessment: Looks good.
- Highest-risk area: None identified.
- Merge recommendation: Safe to merge.

## Blocking

None.

## Important

None.

## Nits

None.

## Test gaps

None.

## EEGLAB notes

None.

Checked for correctness bugs, MATLAB/Octave behavior, EEG structure consistency, GUI/history regressions, changed-behavior tests, and AGENTS.md compliance.
```

Post or update the review comment with:

```bash
gh pr comment "$PR_NUMBER" --edit-last --create-if-none --body-file <file>
```

## Inline Comments

Use inline comments only for validated issues that are best attached to a
specific changed line. Do not duplicate summary findings inline.

Inline comments must:

- State severity.
- Explain the concrete bug or risk.
- Suggest a fix.
- Include a committable suggestion only when the suggestion fully fixes the
  issue and is small.

Never credit yourself or AI tools in review comments.
