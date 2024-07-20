# Contributing to EEGLAB

There are several ways in which you can contribute to the ongoing development of EEGLAB:

## Bug reporting

1. Please search on both the [EEGLAB discussion list archive](https://sccn.ucsd.edu/pipermail/eeglablist/)
   (to search the list Google "eeglablist your question")
   and the [GitHub issue list](https://github.com/sccn/eeglab/issues)
   to see if anybody else has lodged a similar observation.

3. How confident are you that the behavior you have observed is in fact a
   genuine bug, and not a misunderstanding?

   -  *Confident*: Please [open a new GitHub issue](https://github.com/sccn/eeglab/issues/new);
      select the "bug report" issue template to get started.

   -  *Not so confident*: That's fine! Consider instead creating a new topic
      on the [EEGLAB discussion list](https://eeglab.org/others/EEGLAB_mailing_lists.html);
      others can then comment on your observation and determine the
      appropriate level of escalation.

## Requesting a new feature

Please search the [GitHub issue list](https://github.com/sccn/eeglab/issues)
to see if anybody else has made a comparable request:

   -  If a corresponding issue already exists, please add a comment to that
      issue to escalate the request. Additionally, describe any
      aspect of that feature not yet described in the existing issue.

   -  If no such listing exists, then you are welcome to create a [new
      issue](https://github.com/sccn/eeglab/issues/new) outlining the
      request. Be sure to select the "feature request" option to get started
      with writing the issue.

## Asking questions

General questions regarding EEGLAB download, usage, or any other
aspect that is not specific to the EEGLAB *code*, should be directed to
the [EEGLAB discussion list](https://eeglab.org/others/EEGLAB_mailing_lists.html). Also check
the [online documentation](https://eeglab.org/).

## Making direct contributions

Thanks for your interest in making direct contributions to EEGLAB!
We are excited to expand the breadth of researchers involved in improving
and expanding this software, and to ensure that all who make such
contributions receive appropriate acknowledgment through Git.

The instructions below give a short overview of how to go about generating a
proposed change to EEGLAB. A more detailed tutorial on using Git and contributing
to the code (or website) is presented as [online tutorial](https://eeglab.org/tutorials/contribute/)
on the EEGLAB website.

1. You will need to create a *fork* of the [EEGLAB repository](https://github.com/sccn/eeglab)
   into your GitHub account, where unlike the main EEGLAB repository,
   you will have full write access to make the requisite changes.

2. Create a Git branch that is named appropriately according to the
   modifications that are being made. The existing code branch on which
   this new derived branch should be based depends on the nature of the
   proposed change (described later below).

3. Generate one or more Git commits that apply your proposed changes to
   the repository:

   -  Individual commits should ideally have a clear singular purpose,
      and not incorporate multiple unrelated changes. If your proposed
      changes involve multiple disparate components, consider breaking
      those changes up into individual commits.

      Conversely, if multiple code changes are logically grouped with /
      linked to one another, these should ideally be integrated into a
      single commit.

   -  Commits should contain an appropriate message that adequately
      describes the change encapsulated within.

      If the change demands a longer description, then the commit message
      should be broken into a synopsis (less than 80 characters) and message
      body, separated by two newline characters (as this enables GitHub to
      parse them appropriately).

      This can be achieved at the command-line as follows:

      `$ git commit -m $'Commit synopsis\n\nHere is a much longer and wordier description of my proposed changes that doesn\'t fit into 80 characters.\nI can even spread the message body across multiple lines.'`

      (Note also the escape character "`\`" necessary for including an
      apostrophe in the message text)

   -  Where relevant, commit messages can also contain references to
      GitHub issues or pull requests (type the "`#`" character followed
      by the issue / PR number), and/or other individual commits (copy
      and paste the first 8-10 characters of the commit hash).

   -  If multiple persons have contributed to the proposed changes, it is
      possible to modify individual Git commits to have [multiple
      authors](https://help.github.com/en/articles/creating-a-commit-with-multiple-authors),
      to ensure that all contributors receive appropriate acknowledgement.

   As a general rule: Git commits and commit messages should be constructed
   in such a way that, at some time in the future, when one is navigating
   through the contribution history, the evolution of the code is as clear
   as possible.

4. Identify the appropriate classification of the change that you propose
   to make, and read the relevant instructions there:

   -  "[**Fix**](#fix)": If the current code behaviour is
      *clearly incorrect*.

   -  "[**Enhancement**](#enhancement)": If the proposed change improves the *performance* or extends the *functionality* of a particular command or process.

   -  "[**Documentation**](#documentation)": If you made some changes to the function description in the help section of the function.

5. Check that your modified code does not prevent EEGLAB from
   passing existing tests, wherever possible (all test files are in the EEGLAB test directory).

6. For code contributions, if possible, a unit test or reproducibility
   test should be added. This not only demonstrates the behavior of the
   proposed code, but will also preclude future regression of the behavior
   of that code.

1. Once completed, a Pull Request should be generated that merges the
   corresponding branch in your forked version of the EEGLAB repository
   into the appropriate target branch of the main EEGLAB repository
   itself:

   -  If your intended changes are complete, and you consider them ready
      to be reviewed by an EEGLAB team member and merged imminently,
      then create a standard Pull Request.

   -  If your changes are ongoing, and you are seeking feedback from the
      EEGLAB developers before completing them, then also create a standard pull
      request. When you modify your code, the pull request will automatically
      be updated.

#### Fix

1. If there does not already exist a [GitHub issue](https://github.com/sccn/eeglab/issues)
   describing the bug, consider reporting the bug as a standalone issue
   prior to progressing further; that way developers can confirm the issue,
   and possibly provide guidance if you intend to resolve the issue yourself.
   Later, the Pull Request incorporating the necessary changes should then
   reference the listed issue (simply add somewhere in the description
   section the "`#`" character followed by the issue number).

2. Bug fixes are merged directly to `master`; as such, modifications to the
   code should be made in a branch that is derived from `master`, and the
   corresponding Pull Request should select `master` as the target branch
   for code merging.

3. A unit test or reproducibility test should ideally be added: such a
   test should fail when executed using the current `master` code, but pass
   when executed with the proposed changes.

#### Enhancement

1. New features, as well as any code changes that extend the functionality of
   EEGLAB, are merged to the `develop` branch, which contains
   all resolved changes since the most recent tag update. As such, any
   proposed changes that fall under this classification should be made
   in a branch that is based off of the `develop` branch, and the corresponding
   Pull Request should select `develop` as the target branch for code merging.

#### Documentation

If you want to contribute to the documentation on the EEGLAB website, please refer to the website's [Github repository](https://github.com/sccn/sccn.github.io).

#### Coding conventions

Please follow the [MATLAB style guidelines](https://www.mathworks.com/matlabcentral/fileexchange/46056-matlab-style-guidelines-2-0) for guidelines on contributing to the EEGLAB code base.

Non-exhaustive summary of coding guidelines:

* 2 space indents; indent using spaces and not tabs
* No spaces between function name and opening parenthesis of argument list
* One space after the comma that separates function arguments
* Vertically align code with spaces in case it improves readability

#### References

This document is based on the excellent CONTRIBUTING.md document from the [MRTRIX repository](https://github.com/MRtrix3/mrtrix3/), and adjusted accordingly. 
