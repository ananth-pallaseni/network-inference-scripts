using NetworkInference
using Distributions
using SmoothingSplines
using Discretizers
using PyPlot

#################### Functions for getting distributions ####################

# Read the data in filpath into a matrix
function get_data_as_mat(filepath::String; delim::Union{Char,Bool} = false)
    f = open(filepath);
    if delim == false
        mat = readdlm(f, skipstart=1);
    else
        mat = readdlm(f, delim, skipstart=1);
    end
    close(f);
    return mat
end


# Creates gene objects from data by first shuffling the order of the measurements for each gene
function extract_genes_shuffle(datamat::Array{Any, 2}; discretizer = "bayesian_blocks",
    estimator = "maximum_likelihood", number_of_bins = 10)
    num_genes = size(datamat, 1)
    shuffGenes = Array{Gene}(num_genes)
    for i in 1:num_genes
        # Shuffle data values
        shuff = shuffle(datamat[i, 2:end]); # First value in each line is always the id, dont shuffle that

        # Append the name value to the front
        shuff = cat(1, datamat[i, 1], shuff);

        # Convert to gene
        shuffGenes[i] = Gene(shuff, discretizer, estimator, number_of_bins);
    end
    return shuffGenes
end

# Polymorphic version of above function
function extract_genes_shuffle(filepath::String; delim::Union{Char,Bool} = false, discretizer = "bayesian_blocks",
    estimator = "maximum_likelihood", number_of_bins = 10)
    mat = get_data_as_mat(filepath, delim=delim)
    return extract_genes_shuffle(mat, discretizer=discretizer, estimator=estimator, number_of_bins=number_of_bins)
end

# Fits a spline to the input values and returns it as a probability distribution
function get_spline_distr(mi_values::Array{Float64, 1}; λ=0.000011)
    # Create histogram from raw data
    gene_hist = PyPlot.plt[:hist](mi_values, normed=true, bins=100);
    plt[:close]()
    hist_counts = gene_hist[1];

    # Get midpoints for each bar
    hist_dx_over_2 = (gene_hist[2][2] - gene_hist[2][1])/2;
    hist_mids = gene_hist[2][1:end-1] + hist_dx_over_2;

    # Fit spline to histogram points
    fhat_distr = fit(SmoothingSpline, hist_mids, hist_counts, λ);

    return fhat_distr
end

# Takes in a path to gene data and outputs a function mapping from mi to null distr value
function get_f0_distr(filepath; delim::Union{Char,Bool} = false, discretizer = "bayesian_blocks", method=MINetworkInference())
    shuff_genes = extract_genes_shuffle(filepath, delim=delim, discretizer=discretizer)
    na_shuff = NetworkAnalysis(method, shuff_genes)
    mi_shuff = [e.confidence for e in na_shuff.edges]

    # If using CLR, then add epsilon to all values to compensate for the 0 values
    if method isa CLRNetworkInference
        mi_shuff += eps()
    end
    
    f0_distr = fit(Gamma, mi_shuff)
    f0(x) = pdf(f0_distr, x)
    return f0
end

# TAkes in a path to gene data and outputs a function mapping from mi to data distr value
function get_fhat_distr(filepath; delim::Union{Char,Bool} = false, discretizer = "bayesian_blocks", method=MINetworkInference(), λ=0.000011)
    genes = get_genes(filepath, delim=delim, discretizer=discretizer)
    na = NetworkAnalysis(method, genes)
    mi = [e.confidence for e in na.edges]
    fhat_distr = get_spline_distr(mi, λ=λ)
    fhat(x) = max(0, SmoothingSplines.predict(fhat_distr, x))
    return fhat
end

function get_f0_and_fhat_distrs(filepath; delim::Union{Char,Bool} = false, discretizer = "bayesian_blocks", method=MINetworkInference(), λ::Float64=0.000011)
    f0 = get_f0_distr(filepath, delim=delim, discretizer=discretizer, method=method)
    fhat = get_fhat_distr(filepath, delim=delim, discretizer=discretizer, method=method, λ=λ)
    return f0, fhat
end
