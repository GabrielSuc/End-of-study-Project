function offset = compute_offset_for_memorized_buffer(stage, bestPerm, delay_coefficient)
    
    offset = 0;
    for i = 1:(stage - 1)
        offset = offset + delay_coefficient(i,1)*(bestPerm(1,i)*bestPerm(i + length(bestPerm)/2) - 1);
    end
    
end    