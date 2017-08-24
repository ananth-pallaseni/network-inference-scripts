# Grab distribution functions
include("p1_distributions.jl")

#################### Functions for creating networks ####################

function getKeyFromStrings{T <: AbstractString}(s1::T, s2::T)
    if s1 > s2
        tup = Tuple([s2, s1])
    else
        tup = Tuple([s1, s2])
    end
    return tup
end

# Take in an Edge object and output a tuple which can act as a key in a dict
function getKeyFromEdge(e::NetworkInference.Edge)
    gene_names = map(q->q.name, e.genes)
    tup = getKeyFromStrings(gene_names...)
    return tup
end


function mi_to_p1(mi, f0, fhat; prior=1, prior_pos=0.2, prior_neg=0.8)
    p0 = prior == 1 ? prior_pos : prior_neg
    fdr = (f0(mi) * p0) / fhat(mi)
    p1 = 1 - fdr

    if (prior == 1)
        println("mi: ", mi, "    p0: ", p0, "   fdr: ", fdr, "  p1: ", p1, "\n")
    end

  return p1
end



# Creates a network with the top edges by p1 score from the given network
#
# Arguments:
# - na is a NetworkAnalysis object
# - f0 is a function that maps from mi to null distribution value
# - fhat is a function that maps from mi to the data distribution value
# - priors is a dictionary of (Tuple{String, String}, 1/0) representing whether an edge exists or not in the prior data
# - top_threshold is a Float64 of what percentage of top edges to use
#
# Returns a NetworkAnalysis object
function create_network_using_priors(na::NetworkAnalysis, f0, fhat, priors, top_threshold; prior_pos=0.2, prior_neg=0.8)
    # Declare edge array
    p1_edges = Array{Edge, 1}(length(na.edges))

    # Go through edges
    for i in 1:length(na.edges)
        edge = na.edges[i]

        # Get prior value
        key = getKeyFromEdge(edge)
        p0 = get(priors, key, 0)

        if (p0 == 1)
            println(key, ": ", p0)
        end

        # Get p1 value
        p1 = mi_to_p1(edge.confidence, f0, fhat, prior=p0, prior_pos=prior_pos, prior_neg=prior_neg)
        new_edge = Edge(na.edges[i].genes, p1)
        p1_edges[i] = new_edge
    end

    # The mi vs p1 graph looks quadratic, so we only consider the p1 values that came from mi in the top half of the mi distribution. The rest are set to -10
    top_half_of_edges_by_mi = p1_edges[floor(Int64, 1:length(p1_edges)/2)] # mi values are sorted is descending order, so the top half is the first half
    bot_half_of_edges_by_mi = p1_edges[floor(Int64, length(p1_edges)/2):end] # grab bot half
    bot_half_of_edges_by_mi = [Edge(e.genes, -10) for e in bot_half_of_edges_by_mi] # set bot half values to -10
    combined_edges = vcat(bot_half_of_edges_by_mi, top_half_of_edges_by_mi)

    # Sort the edges by p1 score
    sorted_edges = sort(combined_edges, by=e->e.confidence, rev=true)
    sorted_edges = sorted_edges[1:end]
    sorted_genes = map(e->e.genes, sorted_edges)
    sorted_genes = collect(Set(Iterators.flatten(sorted_genes)))
    sorted_na = NetworkAnalysis(sorted_genes, sorted_edges)
    return sorted_na
end



# Creates a matrix such that matrix[i, j] is the network for organism i using prior quality j
#
# Arguments:
# - organism_data_filepaths is an array of filepaths to the gene data
# - prior_qual_matrix is an array of Dictionaries. It has n rows and m columns where n is the number of organisms and m is the number of qualities. Each dictionary has ((Gene1, Gene2), 1/0) pairs representing edges and whether or not they were present in the priors for this quality.
# - threshold_value is the top percentage of edges to take from the network
# - prior_pos is what value to use in lieu of a positive prior value
# - prior_neg is what value to use in lieu of a negative prior value
#
# Returns a matrix of NetworkAnalysis objects
function create_network_matrix(
    organism_data_filepaths::Array{String, 1},
    prior_qual_matrix;
    threshold_value::Float64=1.0,
    prior_pos::Float64=0.2,
    prior_neg::Float64=0.8
    )

    num_organisms = length(organism_data_filepaths)
    num_categories = size(prior_qual_matrix, 2)

    # Declatre matrix
    network_matrix = Array{NetworkAnalysis, 2}(num_organisms, num_categories)

    # Go through data filepaths
    for (i, fname) in enumerate(organism_data_filepaths)
        # Get gene data for organism i
        genes = get_genes(fname)

        # Perform network inference
        na = NetworkAnalysis(MINetworkInference(), genes)

        # Grab prior array for this organism
        organism_priors = prior_qual_matrix[i, 1:end]

        # Grab f0 and fhat distributions
        f0, fhat = get_f0_and_fhat_distrs(fname)

        # Go through prior quality values
        for (j, priors) in enumerate(organism_priors)
            println("PRIORS = ", j)

            # Use priors to create network
            thresh_na = create_network_using_priors(na, f0, fhat, priors, threshold_value; prior_pos=prior_pos, prior_neg=prior_neg)

            # Add network to matrix
            network_matrix[i,j] = thresh_na
        end
    end

    return network_matrix
end
