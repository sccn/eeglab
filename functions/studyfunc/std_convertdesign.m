% std_convertdesign - temporary function converting STUDY design legacy
%                     format to new format.

function STUDY = std_convertdesign(STUDY,ALLEEG);

for index = 1:length(STUDY.design)
    design(index).name              = STUDY.design(index).name;
    design(index).variable(1).label = STUDY.design(index).indvar1;
    design(index).variable(2).label = STUDY.design(index).indvar2;
    design(index).variable(1).value = STUDY.design(index).condition;
    design(index).variable(2).value = STUDY.design(index).group;
    design(index).variable(1).pairing = STUDY.design(index).statvar1;
    design(index).variable(2).pairing = STUDY.design(index).statvar2;
    design(index).cases.label = 'subject';
    design(index).cases.value = STUDY.design(index).subject;
    design(index).include     = STUDY.design(index).includevarlist;
    setinfo = STUDY.design(index).setinfo;
    for c = 1:length(setinfo)
        design(index).cell(c).dataset  = setinfo(c).setindex;
        design(index).cell(c).trials   = setinfo(c).trialindices;
        design(index).cell(c).value    = { setinfo(c).condition setinfo(c).group };
        design(index).cell(c).case     = setinfo(c).subject;
        design(index).cell(c).filebase = setinfo(c).filebase;
    end;
end;

STUDY.design = design;
STUDY = std_selectdesign(STUDY, ALLEEG, STUDY.currentdesign); 
