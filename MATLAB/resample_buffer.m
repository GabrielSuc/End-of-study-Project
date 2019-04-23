function [downsampling_struct_next] = resample_buffer(input, L, M, downsampling_struct_next,buffer_nb , mode)

if mode == 0 %downsampling
    
    quotient = floor(length(input)/M);
    
    if buffer_nb == 1
        first_position_cur  = 1;
        last_position_cur = length(input) - (first_position_cur + quotient*M);
        
        if last_position_cur <= 0
            downsampling_struct_next.downsampled_buffer = zeros(size(input,1), quotient);
            downsampling_struct_next.length_downsampled_buffer = quotient;
        else
            downsampling_struct_next.downsampled_buffer = zeros(size(input,1), quotient + 1);
            downsampling_struct_next.length_downsampled_buffer = quotient + 1;
        end
        
    else
        
        last_position_prev = downsampling_struct_next.last_position;
        if last_position_prev >= 0 
            first_position_cur = M - last_position_prev;
            
        else
            first_position_cur = abs(last_position_prev);
        end
        
        last_position_cur = length(input) - (first_position_cur + quotient*M);
        
        if last_position_cur <= 0
            downsampling_struct_next.downsampled_buffer = zeros(size(input,1), quotient);
            downsampling_struct_next.length_downsampled_buffer = quotient;
        else
            downsampling_struct_next.downsampled_buffer = zeros(size(input,1), quotient + 1);
            downsampling_struct_next.length_downsampled_buffer = quotient + 1;
        end
        
    end
    
    for i = 0:(downsampling_struct_next.length_downsampled_buffer - 1)
        downsampling_struct_next.downsampled_buffer(i + 1) = input(first_position_cur + i*M);
    end    
    
    downsampling_struct_next.last_position = last_position_cur;
end    
end