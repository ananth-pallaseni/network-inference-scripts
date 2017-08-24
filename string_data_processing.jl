"""
Notes on the application of priors:
1) Prior information is obtained as a number
2) It is then scaled to between 0.005 and 0.995. It is not scaled to 0-1 to avoid issues with having priors of 0. By keeping the bottom slightly above 0 we preserve the ability of strong observed information to survive very low priors.
3) The scaled values are inverted by subtracting them from 1. This is because the empirical bayes formula requires the "probability" that an edge is NOT present.
4) All edges are given a uniform base prior of 0.005, and if a prior exists for an edge, it is added to this base value. Thus final priors can be between 0.005 and 1.

Note that scaling the priors does not affect the outcome of the inference as only relative values matter. Applying a prior of 20 to one edge and 10 to the other is the same as applying 2 and 1 respectively.

Note that this method of using priors treats any STRING evidence as positive, so even a low STRING score is some evidence for existence and is better than no STRING score.
Thus edges can be improved by priors, but not penalized by them.
""";


using NetworkInference

# Grab network creation functions
include("p1_networks.jl")

push!(LOAD_PATH, "/cluster/home/avp16/Desktop/data/assortativity_measure")
using Assortativity


raw_data_filepath = "/cluster/home/avp16/Desktop/data/string_data_mouse/mouse_protein_network_data_raw.txt"
data_filepath = "/cluster/home/avp16/Desktop/data/string_data_mouse/mouse_protein_network_data.txt"




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

# returns a dict mapping from name -> id
function read_aliases(aliasfp)
  mat = readdlm(aliasfp, '\t', String)[:, 1:2] # only want id, name columns
  d = Dict()
  for i in 1:size(mat,1)
    v, k = mat[i, :]
    k = uppercase(k)
    v = uppercase(v)
    d[k] = v
  end
  return d
end

# returns a dict mapping from id -> name
function read_aliases_rev(aliasfp)
  mat = readdlm(aliasfp, '\t', String)[:, 1:2] # only want id, name columns
  d = Dict()
  for i in 1:size(mat,1)
    k, v = mat[i, :]
    k = uppercase(k)
    v = uppercase(v)
    d[k] = v
  end
  return d
end


function pairToKey(p1, p2)
  if p1 > p2
    return (p1, p2)
  else
    return (p2, p1)
  end
end

function edgeToIdKey(e, d)
  p1 = d[uppercase(e.genes[1].name)]
  p2 = d[uppercase(e.genes[2].name)]
  return pairToKey(p1, p2)
end

############ Convert protein network to dict, mapping from (id1, id2)->edge_strength
ppi_mat = readdlm(raw_data_filepath, skipstart=1)
edge_dict = Dict()
for i in 1:size(ppi_mat, 1)
  k = pairToKey(ppi_mat[i, 1], ppi_mat[i, 2])
  edge_dict[k] = ppi_mat[i, 3]
end

scale_prior(pv, mv) =  (((pv/mv) - 0.5)*0.99)+0.5 # Scales prior between 0.005 and 0.995
scale_and_invert_prior(pv, mv) = 1 - scale_prior(pv, mv) # inverts prior

# Normalize values
max_prior_val = maximum(values(edge_dict))
edge_dict = map(kv->Pair(kv[1], scale_and_invert_prior(kv[2], max_prior_val)), edge_dict) # All edge confidences are now inverted and scaled between [0.505, 1.495]

############ Create name dictionary, mapping from name->id
alias_filepath = "/cluster/home/avp16/Desktop/data/string_data_mouse/mouse_protein_to_gene_aliases_uniprot_only.txt"

# name dict maps from name->id
name_dict = read_aliases(alias_filepath)

# maps from id -> name
id_dict = read_aliases_rev(alias_filepath)

############ Grab single_cell_data
sc_filepath = "/cluster/home/avp16/Desktop/data/kirsten_data/kirsten_sc_data.csv"

getPrior(e) = get(edge_dict, edgeToIdKey(e, name_dict), 0) + 0.005 # Grabs the prior value + 0.005 for an edge, or 0.005 if it doesnt exist

sc_genes = get_genes(sc_filepath, delim=',', discretizer="uniform_width")

na_mi = NetworkAnalysis(MINetworkInference(), sc_genes)
f0mi, fhmi = get_f0_and_fhat_distrs(sc_filepath, delim=',', discretizer="uniform_width", method=MINetworkInference())
fdrs_mi = get_fdrs(na_mi, f0mi, fhmi)
emp_bayes_mi = map(e->Edge(e.genes, 1-e.confidence*getPrior(e)), fdrs_mi)
bayes_na_mi = NetworkAnalysis(sc_genes, emp_bayes_mi)

na_clr = NetworkAnalysis(CLRNetworkInference(), sc_genes)
f0clr, fhclr = get_f0_and_fhat_distrs(sc_filepath, delim=',', discretizer="uniform_width", method=CLRNetworkInference())
fdrs_clr = get_fdrs(na_clr, f0clr, fhclr)
emp_bayes_clr = map(e->Edge(e.genes, 1-e.confidence*getPrior(e)), fdrs_clr)
bayes_na_clr = NetworkAnalysis(sc_genes, emp_bayes_clr)

na_puc = NetworkAnalysis(PUCNetworkInference(), sc_genes)
f0puc, fhpuc = get_f0_and_fhat_distrs(sc_filepath, delim=',', discretizer="uniform_width", method=PUCNetworkInference())
fdrs_puc = get_fdrs(na_puc, f0puc, fhpuc)
emp_bayes_puc = map(e->Edge(e.genes, 1-e.confidence*getPrior(e)), fdrs_puc)
bayes_na_puc = NetworkAnalysis(sc_genes, emp_bayes_puc)

na_pidc = NetworkAnalysis(PIDCNetworkInference(), sc_genes)
f0pidc, fhpidc = get_f0_and_fhat_distrs(sc_filepath, delim=',', discretizer="uniform_width", method=PIDCNetworkInference())
fdrs_pidc = get_fdrs(na_pidc, f0pidc, fhpidc)
emp_bayes_pidc = map(e->Edge(e.genes, 1-e.confidence*getPrior(e)), fdrs_pidc)
bayes_na_pidc = NetworkAnalysis(sc_genes, emp_bayes_pidc)

inf_types = ["MI", "CLR", "PUC", "PIDC"]
asst_label_path = "/cluster/home/avp16/Desktop/data/assortativity_measure/GeneLists_cell_cycle_pluripotency.tsv"
asst_bayes = [assortativity(na, 0.05, asst_label_path)[end] for na in [bayes_na_mi, bayes_na_clr, bayes_na_puc, bayes_na_pidc]]
figure()
title("Assortativity for single cell data using STRING scores as priors\nKirsten Data")
bar(inf_types, asst_bayes)
ylabel("Assortativity")
xlabel("Information Measure")
