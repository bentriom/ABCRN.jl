
repressilator = @network_model begin
    R1: (G1 => G1 + m1, α/(1+p3^n) + α0)
    R2: (G2 => G2 + m2, α/(1+p1^n) + α0)
    R3: (G3 => G3 + m3, α/(1+p2^n) + α0)
    R4: (m1 => m1 + p1, β * m1)
    R5: (m2 => m2 + p2, β * m2)
    R6: (m3 => m3 + p3, β * m3)
    R7: (m1 => 0, m1)
    R8: (m2 => 0, m2)
    R9: (m3 => 0, m3)
end "Repressilator pkg"

set_observed_var!(repressilator, [:m1, :m2, :m3, :p1, :p2, :p3])
set_time_bound!(repressilator, 15.0)

export repressilator

