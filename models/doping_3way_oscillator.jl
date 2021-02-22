
doping_3way_oscillator = @network_model begin
    R1: (A + B => 2B, rA*A*B)
    R2: (B + C => 2C, rB*B*C)
    R3: (C + A => 2A, rC*C*A)
    R4: (DA + C => DA + A, rC*C*DA)
    R5: (DB + A => DB + B, rA*A*DB)
    R6: (DC + B => DC + C, rB*B*DC)
end "Doping3wayOscillator"

set_x0!(doping_3way_oscillator, [:A,:B,:C,:DA,:DB,:DC], [333,333,333,10,10,10])
set_param!(doping_3way_oscillator, [:rA,:rB,:rC], [1.0, 1.0, 1.0])
set_time_bound!(doping_3way_oscillator, 0.2)

