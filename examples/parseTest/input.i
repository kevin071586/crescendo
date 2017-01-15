# ----------------------------------------------------------
# Comments start with a '#' character
# ----------------------------------------------------------

begin block type block_name
end

begin another command block second_block_name
end

begin finite element model test_model
  database = test.g
  youngs modulus = 10.0  #we can have comments after text
  lame parameter 1 = 0.5
end

begin solver parameters set1
end
