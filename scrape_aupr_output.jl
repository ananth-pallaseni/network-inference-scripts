

function get_interesting_lines(line_list)
  species = []
  aupr_lines = []
  auroc_lines = []

  is_auroc = false

  for line in line_list
    if startswith(line, "[1] \"Running simulation:")
      species_name = split(line, "[1] \"Running simulation:")[2]
      push!(species, species_name)
    elseif startswith(line, "[1]")
      if contains(line, "AUROC.svg")
        is_auroc = true
      end
    else
      if is_auroc
        push!(auroc_lines, line)
      else
        push!(aupr_lines, line)
      end
    end
  end

  return aupr_lines, auroc_lines, species
end;

function parse_interesting_lines(aupr_lines, auroc_lines, species)
  categories = aupr_lines[1]
  categories = split(categories) # convert to array

  aupr_data = aupr_lines[2:end] # remove header
  aupr_data = map(v -> split(v)[2:end], aupr_data)
  aupr_data = map(v -> map(vv->parse(Float64, vv[2:end-1]), v), aupr_data)
  aupr_data = hcat(aupr_data...)' # convert to matrix

  auroc_data = auroc_lines[2:end] # remove header
  auroc_data = map(v -> split(v)[2:end], auroc_data)
  auroc_data = map(v -> map(vv->parse(Float64, vv[2:end-1]), v), auroc_data)
  auroc_data = hcat(auroc_data...)' # convert to matrix

  species = map(v->split(v)[1][1:end-1], species) # convert to string array

  return aupr_data, auroc_data, species, categories
end;


function get_aupr_output{T <: AbstractString}(lines::Array{T,1})
  aupr_lines, auroc_lines, species = get_interesting_lines(lines)
  aupr, auroc, species, categories = parse_interesting_lines(aupr_lines, auroc_lines, species)
  return aupr, auroc, species, categories
end;

function get_aupr_output(filepath::String)
  f = open(filepath)
  lines = readlines(f)
  close(f)
  return get_aupr_output(lines)
end;
