
intracellular_viral_infection = @network_model begin
    R1: (N + T => G + T, k1 * T * cn)
    R2: (N + G => T, k2 * G * cn)
    R3: (N + A + T => S + T, k3 * T * cn * ca)
    R4: (T => ∅, k4 * T)
    R5: (S => ∅, k5 * S)
    R6: (G + S => V, k6 * G * S)
end "Intracellular viral infection pkg"

set_x0!(intracellular_viral_infection, [:N, :A, :G, :T, :S, :V], [10000, 10000, 0, 1, 0, 0])
set_param!(intracellular_viral_infection, [:cn, :ca, :k1, :k2, :k3, :k4, :k5, :k6], [1.0, 1.0, 1.0, 0.025, 100.0, 0.25, 0.2, 7.5E-6])
set_time_bound!(intracellular_viral_infection, 200.0)

export intracellular_viral_infection

