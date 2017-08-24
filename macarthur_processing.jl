
# Grab network creation functions
include("p1_networks.jl")

push!(LOAD_PATH, "/cluster/home/avp16/Desktop/data/assortativity_measure")
using Assortativity


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


datapath_sc = "/cluster/home/avp16/Desktop/data/macarthur_data/macarthur_matched_filtered_sc_data.txt"
datapath_pop = "/cluster/home/avp16/Desktop/data/macarthur_data/macarthur_matched_filtered_pop_data.txt"
genes_sc = get_genes(datapath_sc, discretizer="uniform_width")
genes_pop = get_genes(datapath_pop, discretizer="uniform_width")

asst_label_path = ""



inf_types = ["MI", "CLR", "PUC", "PIDC"]

networks_sc = [
NetworkAnalysis(MINetworkInference(), genes_sc),
NetworkAnalysis(CLRNetworkInference(), genes_sc),
NetworkAnalysis(PUCNetworkInference(), genes_sc),
NetworkAnalysis(PIDCNetworkInference(), genes_sc)
]


networks_pop = [
NetworkAnalysis(MINetworkInference(), genes_pop),
NetworkAnalysis(CLRNetworkInference(), genes_pop),
NetworkAnalysis(PUCNetworkInference(), genes_pop),
NetworkAnalysis(PIDCNetworkInference(), genes_pop)
]


####### Assortativity for pure networks:
function assortativity_pure()
    asst_sc = [assortativity(na, 0.05, asst_label_path)[end] for na in networks_sc]
    asst_pop = [assortativity(na, 0.05, asst_label_path)[end] for na in networks_pop]

    figure()
    title("Assortativity using pure information measures for MacArthur Single Cell data")
    bar(inf_types, asst_sc)
    ylabel("Assortativity")
    xlabel("Information Measure")

    figure()
    title("Assortativity using pure information measures for MacArthur Population data")
    bar(inf_types, asst_pop)
    ylabel("Assortativity")
    xlabel("Information Measure")
end

####### Assortativity for bayes networks - using pop as priors for sc:
function assortativity_sc_base_pop_prior()
    distrs = [
        get_f0_and_fhat_distrs(datapath_sc, discretizer="uniform_width", method=MINetworkInference()),
        get_f0_and_fhat_distrs(datapath_sc, discretizer="uniform_width", method=CLRNetworkInference()),
        get_f0_and_fhat_distrs(datapath_sc, discretizer="uniform_width", method=PUCNetworkInference()),
        get_f0_and_fhat_distrs(datapath_sc, discretizer="uniform_width", method=PIDCNetworkInference())
        ]

    function get_priors_from_pop(pop_network)
        prior_data = Dict()
        for e in pop_network.edges
            key = getKeyFromEdge(e)
            prior_data[key] = e.confidence
        end
        return prior_data
    end

    fdrs = map(q -> get_fdrs(q[1], q[2]...), zip(networks_sc, distrs))

    # grab priors and scale them (scaling allows us to invert the priors easily)
    # Need to inevrt because the prior information should actually be the probability (ish) of an edge NOT existing.
    # Normalizing does not change the final result, but makes it cleaner as the equation expects values in [0.0, 1.0]
    # Furthermore scaling to a range less than [0.0, 1.0] (in this case [0.005, 0.995]) is to make sure that there are no explicit 0.0 values that happen when we invert
    # (as we dont really know what the true theoretical maximum of the prior is we approximate it with the sample maximum, and the scaling allows to preserve the small likelihoods that get destroyed when we do this)
    priors = [get_priors_from_pop(na) for na in networks_pop]
    max_prior_val = [maximum(values(plst)) for plst in priors]
    priors = [map(kv->Pair(kv[1], 1-((((kv[2]/mv)-0.5)*0.99)+0.5)), pdict) for (pdict, mv) in zip(priors, max_prior_val)] # scale to between 0.005 and 0.995, then invert

    bayes = [map(e->Edge(e.genes, 1 - (e.confidence * pdict[getKeyFromEdge(e)])), fdrlst) for (fdrlst, pdict) in zip(fdrs, priors)]
    bayes = [sort(blist, by=e->e.confidence, rev=true) for blist in bayes]

    networks_bayes = [NetworkAnalysis(na.genes, edges) for (na, edges) in zip(networks_sc, bayes)]

    asst_bayes = [assortativity(na, 0.05, asst_label_path)[end] for na in networks_bayes]

    figure()
    title("Assortativity for single cell data using population data as priors\nMacarthur Data")
    bar(inf_types, asst_bayes)
    ylabel("Assortativity")
    xlabel("Information Measure")
end

####### Assortativity for bayes networks - using sc as priors for pop:

function assortativity_pop_base_sc_priors()
    distrs = [
        get_f0_and_fhat_distrs(datapath_pop, discretizer="uniform_width", method=MINetworkInference()),
        get_f0_and_fhat_distrs(datapath_pop, discretizer="uniform_width", method=CLRNetworkInference()),
        get_f0_and_fhat_distrs(datapath_pop, discretizer="uniform_width", method=PUCNetworkInference()),
        get_f0_and_fhat_distrs(datapath_pop, discretizer="uniform_width", method=PIDCNetworkInference())
        ]

    function get_priors_from_pop(pop_network)
        prior_data = Dict()
        for e in pop_network.edges
            key = getKeyFromEdge(e)
            prior_data[key] = e.confidence
        end
        return prior_data
    end

    fdrs = map(q -> get_fdrs(q[1], q[2]...), zip(networks_pop, distrs))

    # grab priors and scale them (scaling allows us to invert the priors easily)
    # Need to inevrt because the prior information should actually be the probability (ish) of an edge NOT existing.
    # Normalizing does not change the final result, but makes it cleaner as the equation expects values in [0.0, 1.0]
    # Furthermore scaling to a range less than [0.0, 1.0] (in this case [0.005, 0.995]) is to make sure that there are no explicit 0.0 values that happen when we invert
    # (as we dont really know what the true theoretical maximum of the prior is we approximate it with the sample maximum, and the scaling allows to preserve the small likelihoods that get destroyed when we do this)
    priors = [get_priors_from_pop(na) for na in networks_sc]
    max_prior_val = [maximum(values(plst)) for plst in priors]
    priors = [map(kv->Pair(kv[1], 1-((((kv[2]/mv)-0.5)*0.99)+0.5)), pdict) for (pdict, mv) in zip(priors, max_prior_val)] # scale to between 0.005 and 0.995, then invert

    bayes = [map(e->Edge(e.genes, 1 - (e.confidence * pdict[getKeyFromEdge(e)])), fdrlst) for (fdrlst, pdict) in zip(fdrs, priors)]
    bayes = [sort(blist, by=e->e.confidence, rev=true) for blist in bayes]

    networks_bayes = [NetworkAnalysis(na.genes, edges) for (na, edges) in zip(networks_pop, bayes)]

    asst_bayes = [assortativity(na, 0.05, asst_label_path)[end] for na in networks_bayes]

    figure()
    title("Assortativity for population data using single cell data as priors\nMacarthur Data")
    bar(inf_types, asst_bayes)
    ylabel("Assortativity")
    xlabel("Information Measure")
end


assortativity_pure()
assortativity_sc_base_pop_priors()
assortativity_pop_base_sc_priors()
