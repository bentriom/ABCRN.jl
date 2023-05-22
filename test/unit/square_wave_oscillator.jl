
using ABCRN
load_model("square_wave_oscillator")
swo = square_wave_oscillator

swo_reg1 = deepcopy(swo)
set_param!(swo_reg1, [:p1, :p2, :p3], [1.0, 0.0, 0.0])

swo_reg2 = deepcopy(swo)
set_param!(swo_reg2, [:p1, :p2, :p3], [0.0, 1.0, 0.0])

swo_reg3 = deepcopy(swo)
set_param!(swo_reg3, [:p1, :p2, :p3], [0.0, 0.0, 1.0])

swo_irreg1 = deepcopy(swo)
set_param!(swo_irreg1, [:p1, :p2, :p3], [1.0, 1.0, 0.0])
set_param!(swo_irreg1, [:w1, :w2, :w3], [1.0, 1.0, 0.0])

swo_irreg2 = deepcopy(swo)
set_param!(swo_irreg2, [:p1, :p2, :p3], [0.0, 1.0, 1.0])
set_param!(swo_irreg2, [:w1, :w2, :w3], [0.0, 1.0, 1.0])

swo_irreg3 = deepcopy(swo)
set_param!(swo_irreg3, [:p1, :p2, :p3], [1.0, 0.0, 1.0])
set_param!(swo_irreg3, [:w1, :w2, :w3], [1.0, 0.0, 1.0])

swo_irreg4 = deepcopy(swo)
set_param!(swo_irreg4, [:p1, :p2, :p3], [1.0, 1.0, 1.0])
set_param!(swo_irreg4, [:w1, :w2, :w3], [1.0, 1.0, 1.0])

return true

