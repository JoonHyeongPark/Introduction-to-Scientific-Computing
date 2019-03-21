function gaussian_jordan(square_dimension, matrix, b_vector)
  
  matrix
  b_vector
  augmented_matrix = [matrix, b_vector]
  
  for i = 1 : square_dimension, # initialize row pointer and row max.
    row_pointer(i) = i;
    row_max(i) = max(abs(matrix(i,:)));
    
    if(row_max(i) == 0), # if for some row max = 0, then there is no unique solution. -> shut down
      printf("There is no unique solution.\n");
      return;
    end;
    
  end;
  
  row_max
  
  for i = 1 : square_dimension, # elimination process to make diagonal matrix.
    
    M = -1;
    pointer = -1;
    for j = i : square_dimension, 
      if abs(augmented_matrix(row_pointer(j), i)) / row_max(row_pointer(j)) > M, # scale pivotting, find optimal pivot.
        M = abs(augmented_matrix(row_pointer(j), i)) / row_max(row_pointer(j));
        pointer = j;
      end;
    end;
    
    if augmented_matrix(row_pointer(pointer), i) == 0, # if for some pivot = 0, then there is no unique solution. -> shut down
      printf("There is no unique solution.\n");
      return;
    end;
    
    if row_pointer(pointer) != row_pointer(i), # interchange row by scaled pivotting.
      temp = row_pointer(i);
      row_pointer(i) = row_pointer(pointer);
      row_pointer(pointer) = temp;
    end;
    
    for j = 1 : square_dimension, # eliminate column except diagonal pivot.
      
      if i == j, continue; end;
      
      if augmented_matrix(row_pointer(j), i) != 0, # if target column's value is already 0 then there is no need to excute operation.
        
        pivot = augmented_matrix(row_pointer(i), i);
        
        multiplier = augmented_matrix(row_pointer(j), i) / pivot;
        augmented_matrix(row_pointer(j), i) = 0; # A[NROW[j]][i] = A[NROW[j]][i] - A[NROW[j]][i] / A[NROW[i]][i] * A[NROW[i]][i] is obviously zero.
        augmented_matrix(row_pointer(j), i + 1 : square_dimension + 1) -= multiplier * augmented_matrix(row_pointer(i), i + 1 : square_dimension + 1); # to eliminate i-th column, excute row operation.
        
      end;
      
    end;
    
    printf("%dth column has been eliminated.\nAugmented Matrix = \n", i); # i-th column except pivot is eliminated.
    disp(augmented_matrix(row_pointer,:));
    printf("\n");
    
    # if you want more accuracy, use complete pivotting not partial pivotting.
    # just calculate row max again in here.
    # if you add this code below, you can get more accuracy but less speed using complete pivotting.
    
    # for j = i + 1 : square_dimension, row_max(row_pointer(j)) = max(abs(augmented_matrix(row_pointer(j), i + 1 : square_dimension))); end;
    
  end;

  printf("Gaussian Jordan Method with Scaled Partial Pivotting has been finished.\n");
  
  x_vector = augmented_matrix(row_pointer, square_dimension + 1) ./ diag(augmented_matrix(row_pointer,:)); # our answer, there is no need to excute backward subtraction.
  printf("\nX vector has been calculated.\n");
  printf("\nX : \n");
  disp(x_vector);
  
end;