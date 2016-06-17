function [block] = block_dataterm_sublabel(row, col, nx, ny, L, left, right)
    
    data = { nx, ny, L, left, right };
    block = { { 'dataterm_sublabel', row, col, data }, {nx*ny*(L-1), ...
                        nx*ny*(L-1) } }; 
    
end
