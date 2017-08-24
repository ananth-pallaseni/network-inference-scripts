# Creates a folder heirarchy with the following structure:
# base_folder
#     algorithm_folder
#         orgransim_i
#             gold_standard.txt
#             sub_folder_j
#               network.txt
#             ...
#         ...
#
# Input should be such that:
# - organsim_folder_names[i] matches gold_standards_path_array[i]
# - network_matrix[i][j] is the NetworkAnalysis object that should be inside
#   organsim_folder_names[i] and sub_folder_names[j]
function create_and_fill_directory_structure(
  base_folder_name::String,
  algorithm_folder_name::String,
  organism_folder_names::Array{String,1},
  gold_standards_path_array::Array{String,1},
  sub_folder_names::Array{String,1},
  network_matrix::Array{NetworkAnalysis, 2}
  )

  mkdir(base_folder_name) # create outermost directory

  alg_dir_path = base_folder_name * "/" * algorithm_folder_name
  mkdir(alg_dir_path) # create algorithm dir

  # Create and fill organsim dirs
  for i in 1:length(organism_folder_names)
    # Get organism directory name
    org_dir_name = organism_folder_names[i]
    org_dir_path = alg_dir_path * "/" * org_dir_name
    mkdir(org_dir_path) # Create organism directory

    # Copy over golden standard file
    gs_src_path = gold_standards_path_array[i]
    gs_dst_path = org_dir_path * "/" * org_dir_name * "_goldstandard.tsv"
    cp(gs_src_path, gs_dst_path)

    # Create and fill sub dirs
    for j in 1:length(sub_folder_names)
      # Create sub dir
      sub_dir_name = sub_folder_names[j]
      sub_dir_path = org_dir_path * "/" * sub_dir_name
      mkdir(sub_dir_path) # Create sub directory

      # Create network file
      network_path = sub_dir_path * "/network.txt"
      network = network_matrix[i, j]
      write_network_file(network_path, network)
    end
  end
end


# Does the same as the above function, but handles instances where network_matrix
# is not a 2 dimensional arrray, but an array of arrays.
function create_and_fill_directory_structure(
  base_folder_name::String,
  algorithm_folder_name::String,
  organism_folder_names::Array{String,1},
  gold_standards_path_array::Array{String,1},
  sub_folder_names::Array{String,1},
  network_matrix::Array{Array{NetworkAnalysis,1}, 1}
  )

  # Restack the array as a 2 dimensional one
  # The call to permutedims is an unfortunate necessity as the transpose function does not work on nonnumeric arrays at the moment.
  network_matrix_2_dim = permutedims(hcat(network_matrix...), (2,1))

  create_and_fill_directory_structure(
    base_folder_name,
    algorithm_folder_name,
    organism_folder_names,
    gold_standards_path_array,
    sub_folder_names,
    network_matrix,
    )

end



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
