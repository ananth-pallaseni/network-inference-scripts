# Grab network creation functions
include("p1_networks.jl")

# Grab functions for setting up and dealing with aupr
include("aupr_functions.jl")


# Compute false discovery rate = (f0/fhat) values from the given network
#
# Arguments:
# - na : NetworkAnalysis Object
# - f0 : A function that represents the pdf of the null distribution
# - fh : A function that represents the pdf of the data distribution
#
# Returns:
# A list of edges with fdr values instead of confidence values (in order of initial confidence value)
function get_fdrs(na, f0, fh)
    na_genes = [e.genes for e in na.edges]

    mis = [e.confidence for e in na.edges]

    fdrs = [f0(x)/fh(x) for x in mis]

    # Set the bottom half of the fdr values to 1 as these are the lowest mi values
    mididx = floor(Int64, length(mis)/2)
    fdrs[mididx:end] = 1.0

    # set all values greater than 1 to 1
    fdrs[fdrs .> 1.0] = 1.0

    # Edge list
    fdr_edges = [Edge(tup[1], tup[2]) for tup in zip(na_genes, fdrs)]
end



# Returns a dictionary of (gene1, gene2) -> 1/0
function parse_prior_file(filepath::String)
    f = open(filepath)
    lines = readlines(f)
    close(f)

    lines = map(split, lines)

    d = Dict();

    for (g1, g2, present) in lines
        key = getKeyFromStrings(g1, g2)
        val = parse(Int64, present)

        # The prior file is directed, so it has (g1, g2) as well as (g2, g1) and they may have different values. Our method is undirected, so if an edge exists between either pair, we take it to mean an undirected edge exists between them.
        potentially_existing_value = get(d, key, 0) # Get the value that already exists for this pair, or 0 if it doesnt
        d[key] = max(val, potentially_existing_value) # Take the max of the current and existing values.
    end

    return d
end



IS_EDGE = 0.2
NOT_EDGE = 0.8

############################################################
function apply_priors(fdr_edges, org_name, is_edge, not_edge)
    # 0.1 priors
    pr_01 = parse_prior_file("single_cell_data/priors/" * org_name * "_0.1.tsv")
    fdr_01_edges = [pr_01[getKeyFromEdge(e)] == 1 ? Edge(e.genes, e.confidence * is_edge) : Edge(e.genes, e.confidence * not_edge) for e in fdr_edges]
    p1_01_edges = [Edge(e.genes, 1-e.confidence) for e in fdr_01_edges]


    # 0.5 priors
    pr_05 = parse_prior_file("single_cell_data/priors/" * org_name * "_0.5.tsv")
    fdr_05_edges = [pr_05[getKeyFromEdge(e)] == 1 ? Edge(e.genes, e.confidence * is_edge) : Edge(e.genes, e.confidence * not_edge) for e in fdr_edges]
    p1_05_edges = [Edge(e.genes, 1-e.confidence) for e in fdr_05_edges]


    # 0.9 priors
    pr_09 = parse_prior_file("single_cell_data/priors/" * org_name * "_0.9.tsv")
    fdr_09_edges = [pr_09[getKeyFromEdge(e)] == 1 ? Edge(e.genes, e.confidence * is_edge) : Edge(e.genes, e.confidence * not_edge) for e in fdr_edges]
    p1_09_edges = [Edge(e.genes, 1-e.confidence) for e in fdr_09_edges]


    # 1.0 priors
    pr_10 = parse_prior_file("single_cell_data/priors/" * org_name * ".tsv")
    fdr_10_edges = [pr_10[getKeyFromEdge(e)] == 1 ? Edge(e.genes, e.confidence * is_edge) : Edge(e.genes, e.confidence * not_edge) for e in fdr_edges]
    p1_10_edges = [Edge(e.genes, 1-e.confidence) for e in fdr_10_edges]

    return p1_01_edges, p1_05_edges, p1_09_edges, p1_10_edges
end
############################################################


function remove_previous_results(org_name)
    try
        rm("test_aupr_input/bayesian_blocks_maximum_likelihood/" * org_name * "/no_priors/network.txt")

        rm("test_aupr_input/bayesian_blocks_maximum_likelihood/" * org_name * "/only_fdr/network.txt")

        rm("test_aupr_input/bayesian_blocks_maximum_likelihood/" * org_name * "/prior_0_1/network.txt")

        rm("test_aupr_input/bayesian_blocks_maximum_likelihood/" * org_name * "/prior_0_5/network.txt")

        rm("test_aupr_input/bayesian_blocks_maximum_likelihood/" * org_name * "/prior_0_9/network.txt")

        rm("test_aupr_input/bayesian_blocks_maximum_likelihood/" * org_name * "/prior_1_0/network.txt")
    catch
        println("No files to remove")
    end
end


function write_all_network_files(org_name, na, p1_no_priors_edges, p1_01_edges, p1_05_edges, p1_09_edges, p1_10_edges)
    write_network_file("test_aupr_input/bayesian_blocks_maximum_likelihood/" * org_name * "/no_priors/network.txt", na)

    write_network_file("test_aupr_input/bayesian_blocks_maximum_likelihood/" * org_name * "/only_fdr/network.txt", NetworkAnalysis(na.genes, p1_no_priors_edges))

    write_network_file("test_aupr_input/bayesian_blocks_maximum_likelihood/" * org_name * "/prior_0_1/network.txt", NetworkAnalysis(na.genes, p1_01_edges))

    write_network_file("test_aupr_input/bayesian_blocks_maximum_likelihood/" * org_name * "/prior_0_5/network.txt", NetworkAnalysis(na.genes, p1_05_edges))

    write_network_file("test_aupr_input/bayesian_blocks_maximum_likelihood/" * org_name * "/prior_0_9/network.txt", NetworkAnalysis(na.genes, p1_09_edges))

    write_network_file("test_aupr_input/bayesian_blocks_maximum_likelihood/" * org_name * "/prior_1_0/network.txt", NetworkAnalysis(na.genes, p1_10_edges))
end


function get_aupr_stats()
    lines = readstring(`Rscript aupr_calcs.R test_aupr_input/`)
    split_lines = split(lines, "\n", keep=false)
    aupr, auroc, species, categories =  get_aupr_output(split_lines)
    return aupr, auroc, species, categories
end

function do_process(org_name, is_edge, not_edge; method=MINetworkInference(), λ=0.000011)

    println("Getting genes for " * org_name)
    genes = get_genes("single_cell_data/testtesttest/"*org_name*".txt")

    println("Performing Network Inference")
    na = NetworkAnalysis(method, genes)

    println("Creating null and data distributions")
    f0, fh = get_f0_and_fhat_distrs("single_cell_data/testtesttest/"*org_name*".txt", method=method, λ=λ)

    println("Applying Priors")
    fdr_edges = get_fdrs(na, f0, fh)

    p1_no_priors_edges = [Edge(e.genes, 1-e.confidence) for e in fdr_edges]

    p1_01_edges, p1_05_edges, p1_09_edges, p1_10_edges = apply_priors(fdr_edges, org_name, is_edge, not_edge)

    println("Writing to aupr input format ")
    remove_previous_results(org_name)

    write_all_network_files(org_name, na, p1_no_priors_edges, p1_01_edges, p1_05_edges, p1_09_edges, p1_10_edges)

    println("Running aupr calcs")
    aupr, auroc, species, categories = get_aupr_stats()

    aupr = aupr[1, 1:end]

    return categories, aupr
end

function share_common(org_name, na, fdr_edges, p1_no_priors_edges, is_edge, not_edge)
    p1_01_edges, p1_05_edges, p1_09_edges, p1_10_edges = apply_priors(fdr_edges, org_name, is_edge, not_edge)

    remove_previous_results(org_name)

    write_all_network_files(org_name, na, p1_no_priors_edges, p1_01_edges, p1_05_edges, p1_09_edges, p1_10_edges)

    aupr, auroc, species, categories = get_aupr_stats()

    aupr = aupr[1, 1:end]

    return categories, aupr
end

function get_auprs_for_prior_weights(org_name, edge_weights_arr;method=MINetworkInference(), λ=0.000011)
    genes = get_genes("single_cell_data/testtesttest/"*org_name*".txt")

    na = NetworkAnalysis(MINetworkInference(), genes)

    f0, fh = get_f0_and_fhat_distrs("single_cell_data/testtesttest/"*org_name*".txt", method=method, λ=λ)

    fdr_edges = get_fdrs(na, f0, fh)

    p1_no_priors_edges = [Edge(e.genes, 1-e.confidence) for e in fdr_edges]

    aupr_arr = Array{Any, 1}(length(edge_weights_arr))
    for i=1:length(edge_weights_arr)
        println("Running iteration " * string(i))
        tup = edge_weights_arr[i]
        aupr_arr[i] = share_common(org_name, na, fdr_edges, p1_no_priors_edges, tup...)
    end
    return aupr_arr
end


###############################################################

org_name = "50_yeast1"
is_edge = 0.2
not_edge = 0.8
method=PUCNetworkInference()
λ=100.

# categories, aupr = do_process(org_name, is_edge, not_edge, method=method, λ=λ)
# bar(categories, aupr)


# # Array of weights in the form (is_edge, not_edge). is_edge is always <= not_edge
# edge_weights_arr = []
# for x=0.:0.2:1.6
#     for y=x+0.2:0.2:1.6
#        push!(edge_weights_arr, (x,y))
#    end
# end
#
# # Returns an array of the form [... ([categories], [auprs]) ... ]
# aupr_arr = @time get_auprs_for_prior_weights(org_name, edge_weights_arr)
#
# cats = aupr_arr[1][1]
# auprs = [a[2] for a in aupr_arr]


function plot_null_and_data_distrs(datapath::String; datatitle::String="Kirsten Data", delim::Union{Char,Bool} = false, discretizer = "bayesian_blocks", λ_mi=0.000011, λ_clr=0.001, λ_puc=100., λ_pidc=0.001)
    mi_f0, mi_fh = get_f0_and_fhat_distrs(datapath, method=MINetworkInference(), λ=λ_mi, delim=delim, discretizer=discretizer)
    clr_f0, clr_fh = get_f0_and_fhat_distrs(datapath, method=CLRNetworkInference(), λ=λ_clr, delim=delim, discretizer=discretizer)
    puc_f0, puc_fh = get_f0_and_fhat_distrs(datapath, method=PUCNetworkInference(), λ=λ_puc, delim=delim, discretizer=discretizer)
    pidc_f0, pidc_fh = get_f0_and_fhat_distrs(datapath, method=PIDCNetworkInference(), λ=λ_pidc, delim=delim, discretizer=discretizer)

    figure()
    suptitle("Null and Data distributions for various information measures.\nData: " * datatitle)
    subplot(221)
    title("MI. lambda=0.00001")
    xlabel("MI Values")
    ylabel("Density")
    plot(0:0.01:0.6, mi_f0.(0:0.01:0.6), color="red", label="Null Distr")
    plot(0:0.01:0.6, mi_fh.(0:0.01:0.6), color="blue", label="Data Distr")
    legend()

    subplot(222)
    title("CLR. lambda=0.001")
    xlabel("CLR Values")
    ylabel("Density")
    plot(0:0.1:9, clr_f0.(0:0.1:9), color="red", label="Null Distr")
    plot(0:0.1:9, clr_fh.(0:0.1:9), color="blue", label="Data Distr")
    legend()

    subplot(223)
    title("PUC. lambda=100.")
    xlabel("PUC Values")
    ylabel("Density")
    plot(1:0.5:100, puc_f0.(1:0.5:100), color="red", label="Null Distr")
    plot(1:0.5:100, puc_fh.(1:0.5:100), color="blue", label="Data Distr")
    legend()

    subplot(224)
    title("PIDC. lambda=0.001")
    xlabel("PIDC Values")
    ylabel("Density")
    plot(0:0.01:2.2, pidc_f0.(0:0.01:2.2), color="red", label="Null Distr")
    plot(0:0.01:2.2, pidc_fh.(0:0.01:2.2), color="blue", label="Data Distr")
    legend()
end

function plot_aupr_vs_prior_quality()
    mi_categories, mi_aupr = do_process(org_name, is_edge, not_edge, method=MINetworkInference(), λ=0.000011)
    clr_categories, clr_aupr = do_process(org_name, is_edge, not_edge, method=CLRNetworkInference(), λ=0.001)
    puc_categories, puc_aupr = do_process(org_name, is_edge, not_edge, method=PUCNetworkInference(), λ=100.)
    pidc_categories, pidc_aupr = do_process(org_name, is_edge, not_edge, method=PIDCNetworkInference(), λ=0.001)

    figure()
    suptitle("AUPR vs prior quality for various measures")

    subplot(221)
    title("MI. lambda=0.00001")
    ylabel("AUPR")
    bar(mi_categories, mi_aupr)
    bar(mi_categories, mi_aupr)
    PyPlot.yticks(0:0.1:1)

    subplot(222)
    title("CLR. lambda=0.001")
    ylabel("AUPR")
    bar(clr_categories, clr_aupr)
    bar(clr_categories, clr_aupr)
    PyPlot.yticks(0:0.1:1)

    subplot(223)
    title("PUC. lambda=100.")
    ylabel("AUPR")
    bar(puc_categories, puc_aupr)
    bar(puc_categories, puc_aupr)
    PyPlot.yticks(0:0.1:1)

    subplot(224)
    title("PIDC. lambda=0.001")
    ylabel("AUPR")
    bar(pidc_categories, pidc_aupr)
    bar(pidc_categories, pidc_aupr)
    PyPlot.yticks(0:0.1:1)
end
