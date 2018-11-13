function array3 = linkArray(array1,array2)
    % Use this function to combine two arrays made of groude delays and
    % related frequencies. The goal is to merge the two frequencies columns
    % by ascending order while keeping the groupe delay values associated to
    % each frequency.
    
    % The arrays should have this shape
    
    % | f1 gd1 |
    % | f2 gd2 |
    % |   .    |
    % |   .    |
    % |   .    |
    % | fn gdn |
    
    % Verify that two array have same dimension
    if size(array1) ~= size(array2)
        error('The two input arrays have different sizes');
    end
    
    % Number of lines in each table
    N1 = size(array1,1);
    N2 = size(array2,1);
    
    % Create a thrid array which will be the final one
    array3 = zeros(N1+N2,2);
    
    % Sart sorting the arrays
    i = 1;
    j = 1;
    k = 1;
    while i < N1 && j < N2
        
    % Check if current element of first 
    % array is smaller than current element 
    % of second array. If yes, store first  
    % array element and increment first array 
    % index. Otherwise do same with second array 

        if array1(i,:) < array2(j,:)
            array3(k,:) = array1(i,:); 
            k = k + 1;
            i = i + 1;
        else
            array3(k,:) = array2(j,:); 
            k = k + 1;
            j = j + 1;
        end
    end
    
    % Store remaining elements of first array
    while i <= N1
        array3(k,:) = array1(i,:);
        k = k + 1;
        i = i + 1;
    end    
        
    % Store remaining elements of second array
    while j <= N2
        array3(k,:) = array2(j,:);
        k = k + 1;
        j = j + 1;
    end 
    

%pensez a checker si pas 2 fois meme ligne
end