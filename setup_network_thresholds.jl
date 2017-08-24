"""
# A script that sets up and populates a directory structure
# to be used with the aupr.R script

# Performs network inference for each of the  given species
# and then produces networks containing the top percentage
# of edges (by confidence) for a given set of cutoffs.
# These thresholded networks are put into the correct folder
# structure to run the R script on.

# COMMAND LINE ARGUMENTS:
# -data sets the path to the folder conainting all the species data
# -outdir sets the path+name of the output directory to create
# -goldstandards sets the path to the directory containing the gold standard files (should have the same names as the corresponding data files)
# -percentages sets an array of percentage cutoffs to use
# -v makes the script print out its progress
""";

using NetworkInference
include("aupr_functions.jl")

############### Constants ###############

# Path to the directory containing expresssion data for the species
path_to_species_data = "/cluster/home/avp16/Desktop/data/single_cell_data/data_to_test"

# Path to directory containing gold standard tsv files. Each gold stamdard file must have the same name as its corresponding expression data file. Eg: species1.tsv should match species1.txt
path_to_gold_standards = "/cluster/home/avp16/Desktop/data/single_cell_data/gold_standards"

# Output folder to create
path_to_output_dir = "aupr_input"

# Top % of edges to take
percentage_thresholds = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35]

# Whether or not to print progress messages
verbose = false

# Indentation of the progress messages
indentation = 0

############### End Constants ###############


############### Functions ###############

# Print out an indented message. Always uses the global indentation level
function print_progress(msg)
  print("    " ^ indentation)
  println(msg)
end

# Go from filename to NetworkAnalysis
function filename_to_network(path_to_file)
  if verbose
    print_progress("Performing inference on " * path_to_file)
  end
  genes = get_genes(path_to_file)
  na = NetworkAnalysis(PIDCNetworkInference(), genes)

  return na
end

# Strips the .txt extension from a file and returns it
get_species_name(fname) = split(fname, ".txt")[1]

# Strips the .tsv extension from a file and returns it
get_gold_standard_name(fname) = split(fname, ".tsv")[1]


# Returns a new NetworkAnalysis object containing the top percentage of edges in the original.
function take_top_percentage_of_edges(na::NetworkAnalysis, top_percentage::Number)
  edges_sorted = sort(na.edges, by=e->e.confidence, rev=true)
  end_index = max(0,  floor(Int64, length(edges_sorted)*top_percentage))
  top_edges = edges_sorted[1:end_index]
  return NetworkAnalysis(na.genes, top_edges)
end


# Input is the path to a directory containing single cell data
# Output is a list of (species_name, NetworkAnalysis) tuples
# Species name is the name of the file minus its extension
function create_networks_from_directory(path_to_directory)
  if verbose
    print_progress("Looking for data in " * path_to_directory)
    global indentation
    indentation += 1
  end

  # create list of files in directory:
  file_list = readdir(path_to_directory)

  append_fname(fname) = path_to_directory * "/" * fname

  ret_lst = map(f->(get_species_name(f), filename_to_network(append_fname(f))), file_list)

  if verbose
    indentation -= 1
  end

  return ret_lst
end

# Creates a directory structure in ouput_dir that can be used in aupr_calcs.R
# Arguments are as follows:
# - network_list is list of (species_name::String, na::NetworkAnalysis) pairs
# - output_dir is a string that is the path to the output directory
# - gold_standard_dir is a string that is the path to the gold standard directory
#
# Returns nothing
function create_threshold_directories(network_list, output_dir, gold_standard_dir, percentages)
  base_folder_name = output_dir

  algorithm_folder_name = "bayesian_blocks_maximum_likelihood"

  organism_folder_names = map(tup->convert(String, tup[1]), network_list)

  sub_folder_names = map(perc->"top_"*string(perc), percentages)

  network_matrix = map(
    na->map(
      perc->take_top_percentage_of_edges(na[2], perc),
      percentages
    ),
    network_list
  )

  gold_standards_path_array = map(tup->gold_standard_dir * "/" * tup[1] * ".tsv", network_list)

  create_and_fill_directory_structure(
    base_folder_name,
    algorithm_folder_name,
    organism_folder_names,
    gold_standards_path_array,
    sub_folder_names,
    network_matrix
  )

  nothing
end

# Performs the network inference and directory setup, given paths to the input, ouput and gold_Standard directories
function do_setup(input_data_path, output_data_path, gold_standard_data_path)
  # Run network inference on each species:
  network_list = create_networks_from_directory(input_data_path)

  # Use network list to create and populate the directories we need
  create_threshold_directories(network_list, output_data_path, gold_standard_data_path, percentage_thresholds)
end


# Checks if flag is a keyword, and if it is, then set the appropriate value
function check_flag_and_set(flag, val)
  if flag == "-data"
    global path_to_species_data
    path_to_species_data = val
  elseif flag == "-outdir"
    global path_to_output_dir
    path_to_output_dir = val
  elseif flag == "-goldstandards"
    global path_to_gold_standards
    path_to_gold_standards = val
  elseif flag == "-percentages"
    global percentage_thresholds
    percentage_thresholds = val
  elseif flag == "-v"
    global verbose
    if val == "true"
      verbose = true
    else
      verbose = false
    end
  end
end


# Goes through command line arguments and sets the appropriate variables
function parse_cmd_args()
  i = 1
  while i < length(ARGS)
    check_flag_and_set(ARGS[i], ARGS[i+1])
    i += 2
  end
end

function main()
  # Grab command line arguments
  parse_cmd_args()

  # Perform inference and directory setup
  do_setup(path_to_species_data, path_to_output_dir, path_to_gold_standards)
end

############### END Functions ###############


# For each species, get data, run inference,
#   For each threshold, collect into network.txt
# Put above data into folders
# main()
