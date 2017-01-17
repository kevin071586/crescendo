#------------------------------------------------------
# Initial cut at defining parameters in an input deck
#------------------------------------------------------

begin finite element model test_model
  youngs modulus = 10.0
  poissons ratio = 0.0
  density = 1.0
  database name = two_hex8_elements.g
end

begin results output output_name
  database name = output.e
end


