function [table_array] = convertTable2Array(table_txt)
%CONVERTTABLE2ARRAY converts a table to an array
%
% DESCRIPTION:
%     convertTable2Array converts a table to an array
%
% USAGE:
%     
%
% INPUTS:
%     table_text     - the table in a '.txt' format
%
% OUTPUTS:
%     table_array  - table in the array format
%
% ABOUT:
%     author        - Ashkan Javaherian
%     date          - 05.08.2020
%     last update   - 12.12.2020
%
%
% This function is part of the r-Wave Toolbox (http://www.r-wave.org)
% Copyright (C) 2021 Ashkan Javaherian
%


table_array = zeros(size(table_txt));

for i = 1: size(table_array, 1)
    for j = 1: size(table_array, 2)
        table_array(i, j) = str2double(cell2mat(table_txt(i, j)));
    end
end


end

