function scalval=extractZeissScaling(bfImage)
% The function extracts ROIs from a Zeiss CZI file opened by Bio-Fomats.
% The CZI image file must be opened by the 'bfImage' function of
% Bio-Formats. This is the input to the extractZeissRois function (bfImage). 
% The ROIs are stored in a structure variable (allRois).
% This structure variable can be used by the drawZeissRois function to
% display the ROIs.
%
% Peter Nagy, email: peter.v.nagy@gmail.com, https://peternagy.webs.com/
% V1.öö
scalval=extractScalingInfomationFromZeiss(bfImage{2});
end
function scalval=extractScalingInfomationFromZeiss(hashTable)
allKeys = arrayfun(@char, hashTable.keySet.toArray, 'UniformOutput', false);
allValues = cellfun(@(x) hashTable.get(x), allKeys, 'UniformOutput', false);
roiDefinitionString = 'Global Scaling|Distance|Value #1';

scal = findSubstringInCellArray(allKeys,roiDefinitionString);
scalval = str2double(allValues{scal});

end
function [indices,strings,positionsInStrings]=findSubstringInCellArray(cellArray,substring)
logicals=strfind(cellArray,substring);
indices=find(cellfun(@(x) ~isempty(x),logicals));
strings=cellArray(indices);
positionsInStrings=cell2mat(logicals(indices));
end