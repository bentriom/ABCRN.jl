
get_module_path() = dirname(dirname(pathof(@__MODULE__)))

function create_results_dir()
    return "./"
end

function cosmos_get_values(name_file::String) 
    output_file = open(name_file)
    dict_values = Dict{String}{Float64}()
    while !eof(output_file)
        line = readline(output_file)
        splitted_line = split(line, ':')
        if (length(splitted_line) > 1) && tryparse(Float64, splitted_line[2]) !== nothing 
            dict_values[splitted_line[1]] = parse(Float64, splitted_line[2])
        end
    end
    close(output_file)
    return dict_values
end

load_plots() = include(get_module_path() * "/core/plots.jl")

