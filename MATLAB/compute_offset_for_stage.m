function offset = compute_offset_for_stage(stage, L, M)

    offset = 0;
    for i = 1:(stage - 1)
        offset = offset + L(1,i)*M(1,i);
    end



end