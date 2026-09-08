function tests = test_pop_topoplot_colorbar
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'functions'), 'IncludingSubfolders', true));
testCase.TestData.EEG = pop_loadset('filename', 'eeglab_data_epochs_ica.set', 'filepath', fullfile(root, 'sample_data'));
end

function setup(testCase)
testCase.TestData.figures = findall(groot, 'Type', 'figure');
testCase.TestData.visibility = get(groot, 'DefaultFigureVisible');
set(groot, 'DefaultFigureVisible', 'off');
figure('Visible', 'off');
end

function teardown(testCase)
delete(setdiff(findall(groot, 'Type', 'figure'), testCase.TestData.figures));
set(groot, 'DefaultFigureVisible', testCase.TestData.visibility);
end

function testPositiveLimits(testCase)
checkScale(testCase, [1 2], 2, false);
end

function testNegativeLimits(testCase)
checkScale(testCase, [-2 -1], 2, false);
end

function testSymmetricLimits(testCase)
checkScale(testCase, [-2 2], 2, true);
end

function testAsymmetricLimits(testCase)
checkScale(testCase, [-1 3], 2, true);
end

function testZeroLowerEndpoint(testCase)
checkScale(testCase, [0 2], 2, false);
end

function testZeroUpperEndpoint(testCase)
checkScale(testCase, [-2 0], 2, false);
end

function testMultiplePositiveMaps(testCase)
checkScale(testCase, [1 2], [1 2], false);
end

function testMultipleSymmetricMaps(testCase)
checkScale(testCase, [-2 2], [1 2], true);
end

function testDefaultLimits(testCase)
pop_topoplot(testCase.TestData.EEG, 0, 2, 'Component', [], 0);
bar = findall(gcf, 'Type', 'axes', 'Tag', 'cbar');
verifyNumElements(testCase, bar, 1);
verifyTrue(testCase, all(diff(get(bar, 'YTick')) > 0));
verifyEqual(testCase, cellstr(get(bar, 'YTickLabel')), {'-'; '0'; '+'});
end

function testZeroComponent(testCase)
EEG = testCase.TestData.EEG;
EEG.icawinv(:, 2) = 0;
pop_topoplot(EEG, 0, 2, 'Zero component', [], 0);
bar = findall(gcf, 'Type', 'axes', 'Tag', 'cbar');
verifyNumElements(testCase, bar, 1);
verifyTrue(testCase, all(diff(get(bar, 'YTick')) > 0));
end

function testUnrelatedAxesUnchanged(testCase)
other = axes('YTick', [10 20 30]);
figure('Visible', 'off');
checkScale(testCase, [-2 2], 2, true);
verifyEqual(testCase, get(other, 'YTick'), [10 20 30]);
end

function checkScale(testCase, limits, components, signed)
pop_topoplot(testCase.TestData.EEG, 0, components, 'Component', [], 0, 'maplimits', limits);
bar = findall(gcf, 'Type', 'axes', 'Tag', 'cbar');
verifyNumElements(testCase, bar, 1);
ticks = get(bar, 'YTick');
verifyTrue(testCase, all(isfinite(ticks)) && all(diff(ticks) > 0));
labels = cellstr(get(bar, 'YTickLabel'));
if signed
  verifyEqual(testCase, labels, {'-'; '0'; '+'});
  range = get(bar, 'YLim');
  mappedZero = limits(1) + (ticks(2)-range(1))/diff(range)*diff(limits);
  verifyEqual(testCase, mappedZero, 0, 'AbsTol', 1e-12);
else
  values = str2double(labels);
  verifyTrue(testCase, all(isfinite(values)));
  verifyEqual(testCase, values([1 end]), limits(:), 'AbsTol', 1e-12);
end
end
