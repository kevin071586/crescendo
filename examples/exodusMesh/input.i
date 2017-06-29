#------------------------------------------------------
# Initial cut at defining parameters in an input deck
#------------------------------------------------------

begin finite element model test_model
  youngs modulus = 10.0
  poissons ratio = 0.0
  density = 1.0
  database name = two_hex8_elements.g
  #database name = hex8_10x10x10.g
  #database name = hex8_10x10x10_thin.g
end

begin results output output_name
  database name = output.e
end

begin eigen solver my_solver
  #number of modes = 10
  #block size = 30 
  #number of blocks = 5
  #maximum restarts = 100
  #target residual = 1.0e-8

  # small problem settings
  number of modes = 10
  block size = 10
  number of blocks = 2
  maximum restarts = 50
  target residual = 1e-8
end

