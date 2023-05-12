function allRois=extractZeissRois(bfImage)
% The function extracts ROIs from a Zeiss CZI file opened by Bio-Fomats.
% The CZI image file must be opened by the 'bfImage' function of
% Bio-Formats. This is the input to the extractZeissRois function (bfImage). 
% The ROIs are stored in a structure variable (allRois).
% This structure variable can be used by the drawZeissRois function to
% display the ROIs.
%
% Peter Nagy, email: peter.v.nagy@gmail.com, https://peternagy.webs.com/
% V1.öö
allRois=extractRoiInfomationFromZeiss(bfImage{2});
function allRois=extractRoiInfomationFromZeiss(hashTable)
allKeys = arrayfun(@char, hashTable.keySet.toArray, 'UniformOutput', false);
allValues = cellfun(@(x) hashTable.get(x), allKeys, 'UniformOutput', false);
% define ROI definitions
roiTypes={'Rectangle','Circle','Ellipse','Polygon','Bezier'};
roiDefinitionStrings.Rectangle{1}='Global Layer|Rectangle|Geometry|Top';
roiDefinitionStrings.Rectangle{2}='Global Layer|Rectangle|Geometry|Left';
roiDefinitionStrings.Rectangle{3}='Global Layer|Rectangle|Geometry|Width';
roiDefinitionStrings.Rectangle{4}='Global Layer|Rectangle|Geometry|Height';
roiDefinitionStrings.Circle{1}='Global Layer|Circle|Geometry|CenterX';
roiDefinitionStrings.Circle{2}='Global Layer|Circle|Geometry|CenterY';
roiDefinitionStrings.Circle{3}='Global Layer|Circle|Geometry|Radius';
roiDefinitionStrings.Ellipse{1}='Global Layer|Ellipse|Geometry|CenterX';
roiDefinitionStrings.Ellipse{2}='Global Layer|Ellipse|Geometry|CenterY';
roiDefinitionStrings.Ellipse{3}='Global Layer|Ellipse|Geometry|RadiusX'; 
roiDefinitionStrings.Ellipse{4}='Global Layer|Ellipse|Geometry|RadiusY'; 
roiDefinitionStrings.Ellipse{5}='Global Layer|Ellipse|Rotation';
roiDefinitionStrings.Polygon{1}='Global Layer|ClosedPolyline|Geometry|Points';
roiDefinitionStrings.Bezier{1}='Global Layer|ClosedBezier|Geometry|Points';
% search for all recognized ROI types
R=cellfun(@(x) roiDefinitionStrings.(x),roiTypes,'UniformOutput',false);
searchFields=cellfun(@(x) x{1},R,'uniformoutput',false);
q=cellfun(@(x) findSubstringInCellArray(allKeys,x),searchFields,'uniformoutput',false);
numOfRois=numel(vertcat(q{:}));
if numOfRois>0
    allRois(numOfRois,1)=struct;
    % search for individual ROIs
    roiCounter=0;
    for i=1:numel(roiTypes)
        for j=1:numel(roiDefinitionStrings.(roiTypes{i}))
            [indices,strings]=findSubstringInCellArray(allKeys,roiDefinitionStrings.(roiTypes{i}){j});
            posOfHash=cellfun(@(x) strfind(x,'#'),strings);
            localCounter=numel(posOfHash);
            parameterName=strfind(roiDefinitionStrings.(roiTypes{i}){j},'|');
            parameterName=roiDefinitionStrings.(roiTypes{i}){j}(parameterName(end)+1:end);
            for k=1:localCounter
                localNumber=str2double(strings{k}(posOfHash(k)+1:end));
                allRois(roiCounter+localNumber).roiType=roiTypes{i};
                switch roiTypes{i}
                    case {'Rectangle','Circle','Ellipse'}
                        allRois(roiCounter+localNumber).roiData.(parameterName)=str2double(allValues{indices(k)});
                    case {'Polygon','Bezier'}
                        currDir=pwd; % in order to ensure that 'split' does not execute the dipimage version of split
                        cd([matlabroot,filesep,'toolbox',filesep,'matlab',filesep,'strfun']);
                        pointXY=split(allValues{indices(k)},{',',' '});
                        cd(currDir);
                        pointXY=cellfun(@(x) str2double(x),pointXY);
                        pointXY=reshape(pointXY,[2 numel(pointXY)/2])';
                        allRois(roiCounter+localNumber).roiData.(parameterName)=pointXY;
                end
            end
        end
        roiCounter=roiCounter+localCounter;
    end
else
    allRois=[];
end
function [indices,strings,positionsInStrings]=findSubstringInCellArray(cellArray,substring)
logicals=strfind(cellArray,substring);
indices=find(cellfun(@(x) ~isempty(x),logicals));
strings=cellArray(indices);
positionsInStrings=cell2mat(logicals(indices));