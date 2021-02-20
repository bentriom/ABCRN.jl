
repressilator = @network_model begin
    tr1: (G1 => G1 + mRNA1, α/(1+P3^n) + α0)
    tr2: (G2 => G2 + mRNA2, α/(1+P1^n) + α0)
    tr3: (G3 => G3 + mRNA3, α/(1+P2^n) + α0)
    trl1: (mRNA1 => mRNA1 + P1, β * mRNA1)
    trl2: (mRNA2 => mRNA2 + P2, β * mRNA2)
    trl3: (mRNA3 => mRNA3 + P3, β * mRNA3)
    degr1: (mRNA1 => 0, mRNA1)
    degr2: (mRNA2 => 0, mRNA2)
    degr3: (mRNA3 => 0, mRNA3)
    degr4: (P1 => 0, P1)
    degr5: (P2 => 0, P2)
    degr6: (P3 => 0, P3)
end "Repressilator pkg"

set_observed_var!(repressilator, [:mRNA1, :mRNA2, :mRNA3, :P1, :P2, :P3])
set_x0!(repressilator, [:mRNA1, :mRNA2, :mRNA3], fill(0, 3))
set_x0!(repressilator, [:P1, :P2, :P3], [5, 0, 15])
set_param!(repressilator, :n, 2.0)
set_param!(repressilator, [:α, :α0, :β, :n], [400.0, 0.0, 2.0, 2.0])
set_time_bound!(repressilator, 200.0)

