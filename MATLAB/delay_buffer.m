function output = delay_buffer(input, memorized_buffer, length_memorized, delay_coefficient, nbr_of_delay, stage, bestPerm)

    
    output = zeros(size(input,1),length(input) + nbr_of_delay*delay_coefficient(stage,1));

    for i = (length_memorized - nbr_of_delay*delay_coefficient(stage,1) + 1):length_memorized
        output(:,i - (length_memorized - nbr_of_delay*delay_coefficient(stage,1) + 1) + 1) = memorized_buffer(:,i + compute_offset_for_memorized_buffer(stage, bestPerm, delay_coefficient));
    end

    for i = (nbr_of_delay*delay_coefficient(stage, 1) + 1):(length(input) + nbr_of_delay*delay_coefficient(stage, 1))
        output(:,i) = input(:,i - nbr_of_delay*delay_coefficient(stage, 1));
    end

        

        
    


end