function quickTexTable(data,filename,initText,format)

%Turns matrix into latex table format. Skips table header and footer, which
%needs to be written in the Tex file manually. 

%data: matrix of numbers

%filename: string name for output

%initText: optional. adds an extra initial column with title for each row.
%If 'off' does nothing. If any other string then adds that string to every
%row. If a cell of strings size Nrow then adds different title for each row
%Adds an ampersand automatically between title and first column. If want
%ampersand but no text just use ' '

%format: suggest use '%.3f'. The code sets exact 0's to 0, rather than the
%chosen format, but this can be manually edited.




%size of table
Nrows = size(data,1);
Ncols = size(data,2);

%check initText contents
iTisCell = isa(initText,'cell');
iTon = 1;
if ~iTisCell %if text string
    if strcmp(initText,'off') %if off
        iTon = 0; %turn off initial column
    end
end




% open the file
fileID = fopen([filename '.tex'], 'w');

%write the data
for row = 1:Nrows
    %write each column entry followed by &
    if iTon %if initial column on
        if iTisCell %if cell of data, one for each row
            fprintf(fileID,'%s & ',initText{row});
        else %if text string, write same for each row
            fprintf(fileID,'%s & ',initText);
        end
        
    end
    
    
    for col = 1:Ncols
        %write ampersand between cols
        if col > 1
            fprintf(fileID, ' & ');
        end
        %write data
        if abs(data(row,col)) > 1e-10 %this sets 0's to exact 0 in table
            fprintf(fileID, '%s',num2str(data(row,col),format));
        else
            fprintf(fileID, '0');
        end
    end
    
    %add \\ to end of line and move to new line. no \\ added for last line
    if row < Nrows
        fprintf(fileID, ' \\\\ \n');
    end
end

% close the file
fclose(fileID);




