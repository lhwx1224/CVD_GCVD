function myprintmatrix(matrix, fid)
% MYPRINTMATRIX.m prints a given 'matrix' into a text file with the file id
% 'fid'.
% syntax: myprintmatrix(matrix, fid)
% input: matrix - a 2d array
%        fid - file identification variable
%
% Hewenxuan Li 2023 @Cornell

    if fid == -1
        error('Error opening the file for writing!');
    end

    for row = 1:size(matrix, 1)
        for col = 1:size(matrix, 2)
            if isreal(matrix(row, col))
                fprintf(fid, '%f\t', matrix(row, col));
            else
                realPart = round(real(matrix(row, col)), 5);
                imagPart = round(imag(matrix(row, col)), 5);
                if imagPart >= 0
                    fprintf(fid, '%f + %fi\t', realPart, imagPart);
                elseif imagPart < 0
                    fprintf(fid, '%f - %fi\t', realPart, abs(imagPart));
                end
            end
        end
        fprintf(fid, '\n'); % Start a new line after each row
    end
end
