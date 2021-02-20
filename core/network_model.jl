
using MacroTools

function get_multiplicand_and_species(expr::Expr)
    @assert expr.args[1] == :*
    multiplicand = reduce(*, expr.args[2:(end-1)])
    sym_species = expr.args[end]
    return (multiplicand, sym_species)
end
get_multiplicand_and_species(sym::Symbol) = (1, sym)
function get_multiplicand_and_species(expr::Real)
    if expr == 0
        return (0, :∅)
    else
        error("A number can't be a species identifier")
    end
end

#=
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
=#
function get_str_propensity(propensity::Expr, dict_species::Dict, dict_params::Dict)
    operator_expr = propensity.args[1]
    operands_expr = propensity.args[2:end]
    if (operator_expr in [:+, :-]) && length(operands_expr) == 1
        return "($(operator_expr)" * "$(get_str_propensity(operands_expr[1], dict_species, dict_params)))"
    end
    str_propensity = "("
    for op in operands_expr[1:(end-1)]
        str_propensity *= "$(get_str_propensity(op, dict_species, dict_params))" * "$(operator_expr)"
    end
    str_propensity *= "$(get_str_propensity(operands_expr[end], dict_species, dict_params)))"
    return str_propensity
end
function get_str_propensity(propensity::Symbol, dict_species::Dict, dict_params::Dict)
    if haskey(dict_species, propensity)
        return "xn[$(dict_species[propensity])]"
    elseif haskey(dict_params, propensity)
        return "p[$(dict_params[propensity])]"
    else
        error("Error during the parsing of propensity functions: a symbol is neither a parameter or a species.")
    end
end
get_str_propensity(propensity::Real, dict_species::Dict, dict_params::Dict) = "$(propensity)"

function fill_params!(dict_params::Dict{ParameterModel,Int}, l_dim_params::Vector{Int}, 
                      propensity::Expr, list_species::Vector)
    for operand in propensity.args[2:end]
        fill_params!(dict_params, l_dim_params, operand, list_species)
    end
end
function fill_params!(dict_params::Dict{ParameterModel,Int}, l_dim_params::Vector{Int}, 
                      propensity::Symbol, list_species::Vector)
    if !(propensity in list_species) && !haskey(dict_params, propensity)
        l_dim_params[1] += 1
        dict_params[propensity] = l_dim_params[1]
    end
end
fill_params!(dict_params::Dict{ParameterModel,Int}, l_dim_params::Vector{Int}, 
             propensity::Real, list_species::Vector) = nothing

macro network_model(expr_network,expr_name...)
    transitions = Transition[]
    dict_species = Dict{VariableModel,Int}()
    dict_params = Dict{ParameterModel,Int}()
    dim_state = 0
    dim_params = 0
    l_dim_params = [0]
    list_expr_reactions = Any[]
    empty_symbols = [:∅]
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
                        if !haskey(dict_species, sym_species) && !(sym_species in empty_symbols)
                            dim_state += 1
                            dict_species[sym_species] = dim_state
                        end
                    end
                else
                    mult, sym_species = get_multiplicand_and_species(reaction_part)
                    if !haskey(dict_species, sym_species) && !(sym_species in empty_symbols) 
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
    allowed_op_in_propensity = [:*]
    for expr_reaction in list_expr_reactions
        local isreaction = @capture(expr_reaction, TR_: (reactants_ => products_, propensity_))
        fill_params!(dict_params, l_dim_params, propensity, list_species)
        #=
        if typeof(propensity) <: Expr 
            @assert propensity.args[1] in allowed_op_in_propensity "Only product of species/params/constants are allowed in propensity"
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
        =#
    end
    dim_params = l_dim_params[1]
    
    # Creation of names variables
    model_name = isempty(expr_name) ? "Network" : expr_name[1]
    model_name = Symbol(replace(model_name, ' ' => '_') * "Model")
    id = Dates.format(Dates.now(), "YmHMs")
    nbr_reactions = length(list_expr_reactions)
    basename_func = "$(model_name)_$(id)"
    basename_func = replace(basename_func, '-'=>'_')
    
    # Writing of f!
    symbol_func_f! = Symbol("$(basename_func)_f!")
    str_expr_model_f! = "function $(symbol_func_f!)(xnplus1::Vector{Int}, l_t::Vector{Float64}, l_tr::Vector{Transition}, xn::Vector{Int}, tn::Float64, p::Vector{Float64})\n\t"
    # Computation of nu and propensity functions in f!
    str_l_a = "l_a = ("
    str_test_isabsorbing = "@inbounds("
    l_nu = [zeros(Int, dim_state) for i = 1:nbr_reactions]
    for (i, expr_reaction) in enumerate(list_expr_reactions)
        local isreaction = @capture(expr_reaction, TR_: (reactants_ => products_, propensity_))
        # Writing of propensities function
        str_propensity = get_str_propensity(propensity, dict_species, dict_params)
        str_expr_model_f! *= "@inbounds a$(i) = " * str_propensity * "\n\t"
        # Anticipating the write of the function isabsorbing
        str_test_isabsorbing *= str_propensity * "+"
        # Update the nu of the i-th reaction 
        nu = l_nu[i]
        if typeof(reactants) <: Expr && reactants.args[1] == :+ 
            for operand in reactants.args[2:end]
                mult, sym_species = get_multiplicand_and_species(operand)
                if !(sym_species in empty_symbols) 
                    nu[dict_species[sym_species]] -= mult
                end
            end
        else
            mult, sym_species = get_multiplicand_and_species(reactants)
            if !(sym_species in empty_symbols)
                nu[dict_species[sym_species]] -= mult
            end
        end
        if typeof(products) <: Expr && products.args[1] == :+ 
            for operand in products.args[2:end]
                mult, sym_species = get_multiplicand_and_species(operand)
                if !(sym_species in empty_symbols)
                    nu[dict_species[sym_species]] += mult
                end
            end
        else
            mult, sym_species = get_multiplicand_and_species(products)
            if !(sym_species in empty_symbols)
                nu[dict_species[sym_species]] += mult
            end
        end
        str_expr_model_f! *= "nu_$i = $(Tuple(nu))\n\t"
        # Anticipating the line l_a = (..)
        str_l_a *= "a$(i), "
    end
    str_test_isabsorbing = str_test_isabsorbing[1:(end-1)] * ")"
    str_l_a = str_l_a[1:(end-2)] * ")\n\t"
    str_expr_model_f! *= str_l_a
    str_expr_model_f! *= "asum = sum(l_a)\n\t"
    str_expr_model_f! *= "if asum == 0.0\n\t\t"
    str_expr_model_f! *= "copyto!(xnplus1, xn)\n\t\t"
    str_expr_model_f! *= "return nothing\n\t"
    str_expr_model_f! *= "end\n\t"
    # Computation of array of transitions
    str_expr_model_f! *= "l_nu = (" * reduce(*, ["nu_$i, " for i = 1:nbr_reactions])[1:(end-2)] * ")\n\t"
    str_expr_model_f! *= "l_sym_R = $(Tuple(transitions))\n\t"
    # Simulation of the reaction
    str_expr_model_f! *= "u1 = rand()\n\t"
    str_expr_model_f! *= "u2 = rand()\n\t"
    str_expr_model_f! *= "tau = - log(u1) / asum\n\t"
    str_expr_model_f! *= "b_inf = 0.0\n\t" 
    str_expr_model_f! *= "b_sup = a1\n\t" 
    str_expr_model_f! *= "reaction = 0\n\n\t" 
    str_expr_model_f! *= "for i = 1:$(nbr_reactions)\n\t\t"
    str_expr_model_f! *= "if b_inf < asum*u2 < b_sup\n\t\t\t"
    str_expr_model_f! *= "reaction = i\n\t\t\t"
    str_expr_model_f! *= "break\n\t\t"
    str_expr_model_f! *= "end\n\t\t"
    str_expr_model_f! *= "@inbounds b_inf += l_a[i]\n\t\t"
    str_expr_model_f! *= "@inbounds b_sup += l_a[i+1]\n\t"
    str_expr_model_f! *= "end\n\t"
    str_expr_model_f! *= "nu = l_nu[reaction]\n\t"
    str_expr_model_f! *= "for i = 1:$(dim_state)\n\t\t"
    str_expr_model_f! *= "@inbounds xnplus1[i] = xn[i]+nu[i]\n\t"
    str_expr_model_f! *= "end\n\t"
    str_expr_model_f! *= "@inbounds l_t[1] = tn + tau\n\t"
    str_expr_model_f! *= "@inbounds l_tr[1] = l_sym_R[reaction]\n"
    str_expr_model_f! *= "end\n"
   
    # Writing of isabsorbing
    symbol_func_isabsorbing = Symbol("isabsorbing_$(basename_func)")
    str_expr_model_isabsorbing = "$(symbol_func_isabsorbing)(p::Vector{Float64},xn::Vector{Int}) = $(str_test_isabsorbing) === 0.0"
 
    # Creation of code
    expr_model_f! = Meta.parse(str_expr_model_f!)
    expr_model_isabsorbing = Meta.parse(str_expr_model_isabsorbing)

    map_idx_var_model = Dict(value => key for (key, value) in dict_species)
    model_g = [map_idx_var_model[i] for i = 1:length(list_species)]

    return quote
        @everywhere @eval $(MarkovProcesses.generate_code_model_type_def(model_name))
        @everywhere @eval $(MarkovProcesses.generate_code_model_type_constructor(model_name))
        @everywhere @eval $(MarkovProcesses.generate_code_simulation(model_name, symbol_func_f!, symbol_func_isabsorbing))
        @everywhere @eval $expr_model_f!
        @everywhere @eval $expr_model_isabsorbing

        getfield(Main, $(Meta.quot(model_name)))($dim_state, $dim_params, $dict_species, $dict_params, $transitions,
                                                 $(zeros(dim_params)), $(zeros(Int, dim_state)), 0.0, 
                                                 $(Meta.quot(symbol_func_f!)), $(Meta.quot(symbol_func_isabsorbing)); g = $model_g)
    end
end

