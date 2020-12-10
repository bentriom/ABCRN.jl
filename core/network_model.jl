
using MacroTools

function get_multiplicand_and_species(expr::Expr)
    @assert expr.args[1] == :*
    multiplicand = reduce(*, expr.args[2:(end-1)])
    sym_species = expr.args[end]
    return (multiplicand, sym_species)
end
get_multiplicand_and_species(sym::Symbol) = (1, sym)

function get_str_propensity(propensity::Expr, dict_species::Dict, dict_params::Dict)
    str_propensity = ""
    for op in propensity.args[2:end]
        if haskey(dict_species, op)
            str_propensity *= "xn[$(dict_species[op])] * "
        elseif haskey(dict_params, op)
            str_propensity *= "p[$(dict_params[op])] * "
        else
            str_propensity *= "$(op) * "
        end
    end
    return str_propensity[1:(end-2)]
end
function get_str_propensity(propensity::Symbol, dict_species::Dict, dict_params::Dict)
    str_propensity = String(propensity)
    if haskey(dict_species, str_propensity)
        str_propensity = "xn[$(dict_species[str_propensity])]"
    elseif haskey(dict_params, str_propensity)
        str_propensity = "p[$(dict_params[str_propensity])]"
    else
        str_propensity = "$(str_propensity)"
    end
    return str_propensity
end

macro network_model(expr_network,expr_name...)
    model_name = isempty(expr_name) ? "Unnamed macro generated" : expr_name[1]
    transitions = Transition[]
    dict_species = Dict{VariableModel,Int}()
    dict_params = Dict{ParameterModel,Int}()
    dim_state = 0
    dim_params = 0
    list_expr_reactions = Any[]
    # First we detect all of the species
    for expr_reaction in expr_network.args
        local isreaction = @capture(expr_reaction, TR_: (reactants_ => products_, propensity_))
        if isreaction
            push!(list_expr_reactions, expr_reaction)
            push!(transitions, TR)
            # Parsing reactants, products
            for reaction_part in [reactants, products]
                # If there's several species interacting / produced
                if typeof(reaction_part) <: Expr && reaction_part.args[1] == :+ 
                    for operand in reaction_part.args[2:end]
                        mult, sym_species = get_multiplicand_and_species(operand)
                        if !haskey(dict_species, sym_species)
                            dim_state += 1
                            dict_species[sym_species] = dim_state
                        end
                    end
                else
                    mult, sym_species = get_multiplicand_and_species(reaction_part)
                    if !haskey(dict_species, sym_species)
                        dim_state += 1
                        dict_species[sym_species] = dim_state
                    end
                end
            end
        end
        if !isreaction && !(typeof(expr_reaction) <: LineNumberNode)
            error("Error in an expression describing a reaction")
        end
    end
    list_species = [species for species in keys(dict_species)]
    # Then we detect parameters in propensity expressions
    # Parameters are the symbols that are not species (at this point we know all of the involved species)
    for expr_reaction in list_expr_reactions
        local isreaction = @capture(expr_reaction, TR_: (reactants_ => products_, propensity_))
        if typeof(propensity) <: Expr 
            @assert propensity.args[1] == :* "Only product of species/params/constants are allowed in propensity"
            for operand in propensity.args[2:end]
                if typeof(operand) <: Symbol
                    # If it's not a species, it's a parameter
                    if !(operand in list_species) && !haskey(dict_params, operand)
                        dim_params += 1
                        dict_params[operand] = dim_params
                    end
                end
            end
        elseif typeof(propensity) <: Symbol
            if !(propensity in list_species) && !haskey(dict_params, propensity)
                dim_params += 1
                dict_params[propensity] = dim_params
            end
        end
        if !isreaction && !(typeof(expr_reaction) <: LineNumberNode)
            error("Error in an expression describing a reaction")
        end
    end
    # Let's write some lines that creates the function f! (step of a simulation) for this biochemical network
    nbr_rand = rand(1:1000)
    nbr_reactions = length(list_expr_reactions)
    basename_func = "$(replace(model_name, ' '=>'_'))_$(nbr_rand)"
    expr_model_f! = "function $(basename_func)_f!(xnplus1::Vector{Int}, l_t::Vector{Float64}, l_tr::Vector{Transition}, xn::Vector{Int}, tn::Float64, p::Vector{Float64})\n\t"
    # Computation of nu and propensity functions in f!
    str_l_a = "l_a = ("
    str_test_isabsorbing = "@inbounds("
    l_nu = [zeros(Int, dim_state) for i = 1:nbr_reactions]
    for (i, expr_reaction) in enumerate(list_expr_reactions)
        local isreaction = @capture(expr_reaction, TR_: (reactants_ => products_, propensity_))
        # Writing of propensities function
        str_propensity = get_str_propensity(propensity, dict_species, dict_params)
        expr_model_f! *= "@inbounds a$(i) = " * str_propensity * "\n\t"
        # Anticipating the write of the function isabsorbing
        str_test_isabsorbing *= str_propensity * "+"
        # Update the nu of the i-th reaction 
        nu = l_nu[i]
        if typeof(reactants) <: Expr && reactants.args[1] == :+ 
            for operand in reactants.args[2:end]
                mult, sym_species = get_multiplicand_and_species(operand)
                nu[dict_species[sym_species]] -= mult
            end
        else
            mult, sym_species = get_multiplicand_and_species(reactants)
            nu[dict_species[sym_species]] -= mult
        end
        if typeof(products) <: Expr && products.args[1] == :+ 
            for operand in products.args[2:end]
                mult, sym_species = get_multiplicand_and_species(operand)
                nu[dict_species[sym_species]] += mult
            end
        else
            mult, sym_species = get_multiplicand_and_species(products)
            nu[dict_species[sym_species]] += mult
        end
        expr_model_f! *= "nu_$i = $(Tuple(nu))\n\t"
        # Anticipating the line l_a = (..)
        str_l_a *= "a$(i), "
    end
    str_test_isabsorbing = str_test_isabsorbing[1:(end-2)] * ")"
    str_l_a = str_l_a[1:(end-2)] * ")\n\t"
    expr_model_f! *= str_l_a
    expr_model_f! *= "asum = sum(l_a)\n\t"
    expr_model_f! *= "if asum == 0.0\n\t\t"
    expr_model_f! *= "copyto!(xnplus1, xn)\n\t\t"
    expr_model_f! *= "return nothing\n\t"
    expr_model_f! *= "end\n\t"
    # Computation of array of transitions
    expr_model_f! *= "l_nu = (" * reduce(*, ["nu_$i, " for i = 1:nbr_reactions])[1:(end-2)] * ")\n\t"
    expr_model_f! *= "l_sym_R = $(Tuple(transitions))\n\t"
    # Simulation of the reaction
    expr_model_f! *= "u1 = rand()\n\t"
    expr_model_f! *= "u2 = rand()\n\t"
    expr_model_f! *= "tau = - log(u1) / asum\n\t"
    expr_model_f! *= "b_inf = 0.0\n\t" 
    expr_model_f! *= "b_sup = a1\n\t" 
    expr_model_f! *= "reaction = 0\n\n\t" 
    expr_model_f! *= "for i = 1:$(nbr_reactions)\n\t\t"
    expr_model_f! *= "if b_inf < asum*u2 < b_sup\n\t\t\t"
    expr_model_f! *= "reaction = i\n\t\t\t"
    expr_model_f! *= "break\n\t\t"
    expr_model_f! *= "end\n\t\t"
    expr_model_f! *= "@inbounds b_inf += l_a[i]\n\t\t"
    expr_model_f! *= "@inbounds b_sup += l_a[i+1]\n\t"
    expr_model_f! *= "end\n\t"
    expr_model_f! *= "nu = l_nu[reaction]\n\t"
    expr_model_f! *= "for i = 1:$(dim_state)\n\t\t"
    expr_model_f! *= "@inbounds xnplus1[i] = xn[i]+nu[i]\n\t"
    expr_model_f! *= "end\n\t"
    expr_model_f! *= "@inbounds l_t[1] = tn + tau\n\t"
    expr_model_f! *= "@inbounds l_tr[1] = l_sym_R[reaction]\n"
    expr_model_f! *= "end\n"
    expr_model_isabsorbing = "isabsorbing_$(basename_func)(p::Vector{Float64},xn::Vector{Int}) = $(str_test_isabsorbing) === 0.0"
    model_f! = eval(Meta.parse(expr_model_f!))
    model_isabsorbing = eval(Meta.parse(expr_model_isabsorbing))
    return :(ContinuousTimeModel($dim_state, $dim_params, $dict_species, $dict_params, $transitions, 
                                 $(zeros(dim_params)), $(zeros(Int, dim_state)), 0.0, $model_f!, $model_isabsorbing; 
                                 g = $list_species, name=$model_name))
end

