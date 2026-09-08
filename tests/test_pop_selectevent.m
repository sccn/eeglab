function tests = test_pop_selectevent
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'functions'), 'IncludingSubfolders', true));
end

function testRetainsMatchingEpochs(testCase)
EEG = fixture;
selected = pop_selectevent(EEG, 'type', 'target', 'deleteepochs', 'on');
verifyEqual(testCase, selected.trials, 2);
verifyEqual(testCase, selected.data, EEG.data(:, :, [1 3]));
verifyEqual(testCase, numel(selected.event), 4);
verifyEqual(testCase, [selected.event.epoch], [1 1 2 2]);
end

function testEmptySelectionErrorsByDefault(testCase)
verifyError(testCase, @() pop_selectevent(fixture, 'type', 'absent', 'deleteepochs', 'on'), '');
end

function testEmptySelectionAllowed(testCase)
selected = pop_selectevent(fixture, 'type', 'absent', 'deleteepochs', 'on', 'erroronempty', 'off');
verifyEmpty(testCase, selected.data);
verifyEmpty(testCase, selected.event);
end

function testInverseEpochSelection(testCase)
EEG = fixture;
selected = pop_selectevent(EEG, 'type', 'target', 'deleteepochs', 'on', 'invertepochs', 'on');
verifyEqual(testCase, selected.data, EEG.data(:, :, [2 4]));
verifyEqual(testCase, selected.trials, 2);
end

function testDeleteUnselectedEvents(testCase)
selected = pop_selectevent(fixture, 'type', 'target', 'deleteepochs', 'on', 'deleteevents', 'on');
verifyEqual(testCase, {selected.event.type}, {'target' 'target'});
verifyEqual(testCase, [selected.event.epoch], [1 2]);
end

function testKeepEpochsWhenDeletingOnlyEvents(testCase)
EEG = fixture;
selected = pop_selectevent(EEG, 'type', 'target', 'deleteepochs', 'off', 'deleteevents', 'on');
verifyEqual(testCase, selected.data, EEG.data);
verifyEqual(testCase, numel(selected.event), 2);
end

function testExplicitErrorOptionWithNonemptySelection(testCase)
EEG = fixture;
selected = pop_selectevent(EEG, 'type', 'target', 'deleteepochs', 'on', 'erroronempty', 'on');
verifyEqual(testCase, selected.data, EEG.data(:, :, [1 3]));
end

function testDatasetArray(testCase)
EEG = fixture;
selected = pop_selectevent([EEG EEG], 'type', 'target', 'deleteepochs', 'on');
verifyEqual(testCase, numel(selected), 2);
verifyEqual(testCase, [selected.trials], [2 2]);
verifyEqual(testCase, selected(1).data, EEG.data(:, :, [1 3]));
verifyEqual(testCase, selected(2).data, EEG.data(:, :, [1 3]));
end

function EEG = fixture
EEG = eeg_emptyset;
EEG.setname = 'selection regression';
EEG.nbchan = 1;
EEG.srate = 1000;
EEG.pnts = 10;
EEG.trials = 4;
EEG.xmin = 0;
EEG.xmax = 0.009;
EEG.data = reshape(single(1:40), 1, 10, 4);
EEG.event = struct('type', {}, 'latency', {}, 'epoch', {});
for trial = 1:4
  type = 'other';
  if mod(trial, 2), type = 'target'; end
  EEG.event(2*trial-1) = struct('type', type, 'latency', (trial-1)*10+3, 'epoch', trial);
  EEG.event(2*trial) = struct('type', 'distractor', 'latency', (trial-1)*10+7, 'epoch', trial);
end
EEG = eeg_checkset(EEG, 'eventconsistency');
end
