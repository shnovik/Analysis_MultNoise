% Author: Shoshana Novik
% Date: June 28, 2021
% Description: Compiles lists of all up and down state durations from every
%              recording, for re-creating Levenstein et. al. Figure 4a and c 

inputdir = 'C:\\Users\Shoshana\Documents\CSHL Summer\Slow Waves 2019\SlowWaves\';
files = dir(inputdir);

allUpDurs = [];
allDownDurs = [];

for fileIndex = 4:length(files)
    txt = load(fullfile(inputdir, files(fileIndex).name));
    swsStructure = txt.SlowWaves;
    for upIndex = 1:length(swsStructure.ints.UP)
        % Calculate all up durations:
        uDur = swsStructure.ints.UP(upIndex, 2) - swsStructure.ints.UP(upIndex, 1);
        allUpDurs = [allUpDurs; uDur];
    end
    for downIndex = 1:length(swsStructure.ints.DOWN)
        % Calculate all down durations:
        dDur = swsStructure.ints.DOWN(downIndex, 2) - swsStructure.ints.DOWN(downIndex, 1);
        allDownDurs = [allDownDurs; dDur];
    end
end

% Create a dataset in which each row contains two columns: one recorded duration, and in which state it was measured (Up or Down)
% This dataset is then loaded into R
allDurs = [allUpDurs; allDownDurs];
states = [repmat("UP", length(allUpDurs), 1); repmat("DOWN", length(allDownDurs), 1)];
ds = mat2dataset(allDurs);
ds.Properties.VarNames = {'Duration'};
ds.States = nominal(states);
export(ds,'XLSFile','allDurs.xlsx')