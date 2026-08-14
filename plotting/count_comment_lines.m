function n = count_comment_lines(filepath)
%COUNT_COMMENT_LINES Count leading '#' comment lines at the start of a file.
%   n = count_comment_lines(filepath) returns the number of consecutive lines
%   at the beginning of filepath that start with '#'. Stops at the first
%   blank line or non-comment line.
    fid = fopen(filepath, 'r');
    if fid == -1, error('Cannot open file: %s', filepath); end
    n = 0;
    while true
        line = fgetl(fid);
        if ~ischar(line) || isempty(line), break; end
        if line(1) == '#', n = n + 1;
        else, break;
        end
    end
    fclose(fid);
end
