module myapp
using JWAS,DataFrames,CSV,Statistics,Printf
function readapp()
############################################################################
# Read file and save them into a dataframe
############################################################################
file_path = ""
while true
    println("Please enter the file path for your file. \n")
    input = readline()
    file_path = input
    if isfile(file_path)
        println("File found: ", input)
        println("(Press Enter to continue)")
        readline() 
    break 
    else
        println("Unable to locate file. Please try again.\n")
    end
end
println("Press any key to continue")
readline()

# List to store the file data
file_data = []

try
    open(file_path) do file
        for line in eachline(file)
            # Add each line of data to the list
            push!(file_data, line)
        end
    end
    println("File content has been read and stored.")
catch e
    println("An error occurred: ", e)
end

println("Press any key to continue")
readline()

############################################################################
# Read phenotypes into environment
############################################################################

# Find lines starting with #pheno_header, #pheno_separator, and #pheno_missingstrings in pedfile_data
phenofile_line = findfirst(line -> startswith(line, "#phenofile"), file_data)
phenotype_header_line = findfirst(line -> startswith(line, "#phenotype_header"), file_data)
phenotype_separator_line = findfirst(line -> startswith(line, "#phenotype_separator"), file_data)
phenotype_missingstrings_line = findfirst(line -> startswith(line, "#phenotype_missingstrings"), file_data)

# Check if all required lines are found
if phenofile_line !== nothing && phenotype_header_line !== nothing && phenotype_separator_line !== nothing && phenotype_missingstrings_line !== nothing
    phenotype_header_value = strip(split(file_data[phenotype_header_line], ':')[2])
    phenotype_separator_value = strip(split(file_data[phenotype_separator_line], ':')[2])
    phenotype_missingstrings_value = strip(split(file_data[phenotype_missingstrings_line], ':')[2])
    phenofile_value = split(file_data[phenofile_line], ':', limit=2)[2]
    phenofile_value = strip(phenofile_value)

    # Process the values further
    phenofile = string(phenofile_value)
    phenotype_header = parse(Bool, string(phenotype_header_value))
    phenotype_separator = string(phenotype_separator_value)
    phenotype_missingstrings = string(phenotype_missingstrings_value)
    phenotypes = CSV.read(phenofile, DataFrame, delim=phenotype_separator, header=phenotype_header, missingstring=[phenotype_missingstrings])
else
    # Do nothing if any of the required lines is not found
end


############################################################################
# Read pedigrees into environment
############################################################################

# Find lines starting with #pedigree_header, #pedigree_separator, and #pedigree_missingstrings in pedfile_data
pedfile_line = findfirst(line -> startswith(line, "#pedfile"), file_data)
pedigree_header_line = findfirst(line -> startswith(line, "#pedigree_header"), file_data)
pedigree_separator_line = findfirst(line -> startswith(line, "#pedigree_separator"), file_data)
pedigree_missingstrings_line = findfirst(line -> startswith(line, "#pedigree_missingstrings"), file_data)

# Check if all required lines are found
if pedfile_line !== nothing && pedigree_header_line !== nothing &&
    pedigree_separator_line !== nothing && pedigree_missingstrings_line !== nothing
    pedigree_header_value = strip(split(file_data[pedigree_header_line], ':')[2])
    pedigree_separator_value = strip(split(file_data[pedigree_separator_line], ':')[2])
    pedigree_missingstrings_value = strip(split(file_data[pedigree_missingstrings_line], ':')[2])
    pedfile_value = split(file_data[pedfile_line], ':', limit=2)[2]
    pedfile_value = strip(pedfile_value)

    # Process the values further and obtain pedigrees
    pedfile = string(pedfile_value)
    pedigree_header = parse(Bool, string(pedigree_header_value))
    pedigree_separator = string(pedigree_separator_value)
    pedigree_missingstrings = string(pedigree_missingstrings_value)
    pedigree = get_pedigree(pedfile, header=pedigree_header, separator=pedigree_separator, missingstrings=pedigree_missingstrings)
else
    # Do nothing if any of the required lines is not found
end
############################################################################
# Read genotypes into environment
############################################################################
indices = []

for (i, line) in enumerate(file_data)
    if startswith(line, "#genotypes")
        push!(indices, i)
    end
end

g_type = Dict()
genoterm = []


for i in eachindex(indices)
    current_index = indices[i]
    key_str = strip(split(file_data[current_index + 1], ':')[2])
    key = Symbol(key_str)


    push!(genoterm, key)
    
    start_index = current_index + 2
    # Set end_index to the line before the next line starting with #genotypes
    if i < length(indices)
        end_index = indices[i+1] - 1
    else
        # Handle the last case, set end_index to the last line of the file
        end_index = findlast(line -> startswith(line, "#genotype"), file_data)
    end

    # Handle start_index and end_index here, for example, call your process_data_range function
    genofile_line = findfirst(line -> startswith(line, "#genofile"), file_data[start_index:end_index])
    genotype_G_line = findfirst(line -> startswith(line, "#genotype_G"), file_data[start_index:end_index])
    genotype_method_line = findfirst(line -> startswith(line, "#genotype_method"), file_data[start_index:end_index])
    genotype_pi_line = findfirst(line -> startswith(line, "#genotype_pi"), file_data[start_index:end_index])
    genotype_estimatePi_line = findfirst(line -> startswith(line, "#genotype_estimatePi"), file_data[start_index:end_index])
    genotype_estimateVariance_line = findfirst(line -> startswith(line, "#genotype_estimateVariance"), file_data[start_index:end_index])
    genotype_estimateScale_line = findfirst(line -> startswith(line, "#genotype_estimateScale"), file_data[start_index:end_index])
    genotype_separator_line = findfirst(line -> startswith(line, "#genotype_separator"), file_data[start_index:end_index])
    genotype_header_line = findfirst(line -> startswith(line, "#genotype_header"), file_data[start_index:end_index])
    genotype_rowID_line = findfirst(line -> startswith(line, "#genotype_rowID"), file_data[start_index:end_index])
    genotype_center_line = findfirst(line -> startswith(line, "#genotype_center"), file_data[start_index:end_index])
    genotype_G_is_marker_variance_line = findfirst(line -> startswith(line, "#genotype_G_is_marker_variance"), file_data[start_index:end_index])
    genotype_df_line = findfirst(line -> startswith(line, "#genotype_df"), file_data[start_index:end_index])
    genotype_starting_value_line = findfirst(line -> startswith(line, "#genotype_starting_value"), file_data[start_index:end_index])
    genotype_quality_control_line = findfirst(line -> startswith(line, "#genotype_quality_control"), file_data[start_index:end_index])
    genotype_MAF_line = findfirst(line -> startswith(line, "#genotype_MAF"), file_data[start_index:end_index])
    genotype_missing_value_line = findfirst(line -> startswith(line, "#genotype_missing_value"), file_data[start_index:end_index])

    file_data1 = file_data[start_index:end_index]

    if genofile_line !== nothing
        genofile_value = string(strip(split(file_data1[genofile_line], ':', limit=2)[2]))
    else
        println("No line starting with #phenotype found in the file. And they are necessary")
    end
    
    if genotype_G_line !== nothing
        genotype_G_value = strip(split(file_data1[genotype_G_line], ':')[2])
    else
        genotype_G_value = false
    end
    
    if genotype_method_line !== nothing
        genotype_method_value = strip(split(file_data1[genotype_method_line], ':')[2])
    else
        genotype_method_value = "BayesC"
    end
    
    if genotype_pi_line !== nothing
        genotype_pi_value = strip(split(file_data1[genotype_pi_line], ':')[2])
    else
        genotype_pi_value = 0.0
    end
    
    if genotype_estimatePi_line !== nothing
        genotype_estimatePi_value = strip(split(file_data1[genotype_estimatePi_line], ':')[2])
    else
        genotype_estimatePi_value = true
    end
    
    if genotype_estimateVariance_line !== nothing
        genotype_estimateVariance_value = strip(split(file_data1[genotype_estimateVariance_line], ':')[2])
    else
        genotype_estimateVariance_value = true
    end
    
    if genotype_estimateScale_line !== nothing
        genotype_estimateScale_value = strip(split(file_data1[genotype_estimateScale_line], ':')[2])
    else
        genotype_estimateScale_value = false
    end
    
    if genotype_separator_line !== nothing
        genotype_separator_value = collect(strip(split(file_data1[genotype_separator_line], ':')[2]))[1]
    else
        genotype_separator_value = ','
    end
    
    if genotype_header_line !== nothing
        genotype_header_value = strip(split(file_data1[genotype_header_line], ':')[2])
    else
        genotype_header_value = true
    end
    
    if genotype_rowID_line !== nothing
        genotype_rowID_value = strip(split(file_data1[genotype_rowID_line], ':')[2])
    else
        genotype_rowID_value = false
    end
    
    if genotype_center_line !== nothing
        genotype_center_value = strip(split(file_data1[genotype_center_line], ':')[2])
    else
        genotype_center_value = true
    end
    
    if genotype_G_is_marker_variance_line !== nothing
        genotype_G_is_marker_variance_value = strip(split(file_data1[genotype_G_is_marker_variance_line], ':')[2])
    else
        genotype_G_is_marker_variance_value = false
    end
    
    if genotype_df_line !== nothing
        genotype_df_value = strip(split(file_data1[genotype_df_line], ':')[2])
    else
        genotype_df_value = 4.0
    end
    
    if genotype_starting_value_line !== nothing
        genotype_starting_value_value = strip(split(file_data1[genotype_starting_value_line], ':')[2])
    else
        genotype_starting_value_value = false
    end
    
    if genotype_quality_control_line !== nothing
        genotype_quality_control_value = strip(split(file_data1[genotype_quality_control_line], ':')[2])
    else
        genotype_quality_control_value = true
    end
    
    if genotype_MAF_line !== nothing
        genotype_MAF_value = strip(split(file_data1[genotype_MAF_line], ':')[2])
    else
        genotype_MAF_value = 0.01
    end
    
    if genotype_missing_value_line !== nothing
        genotype_missing_value_value = strip(split(file_data1[genotype_missing_value_line], ':')[2])
    else
        genotype_missing_value_value = 9.0
    end

    genotypes = get_genotypes(string(genofile_value),
                        G = genotype_G_value,
                        method = genotype_method_value, 
                        Pi = genotype_pi_value, 
                        estimatePi = genotype_estimatePi_value, 
                        estimateVariance = genotype_estimateVariance_value, 
                        estimateScale = genotype_estimateScale_value, 
                        separator = genotype_separator_value,
                        header = genotype_header_value,
                        rowID = genotype_rowID_value,
                        center = genotype_center_value, 
                        G_is_marker_variance = genotype_G_is_marker_variance_value,
                        df = genotype_df_value, 
                        starting_value = genotype_starting_value_value,
                        quality_control = genotype_quality_control_value, 
                        MAF = genotype_MAF_value, 
                        missing_value = genotype_missing_value_value);

    g_type[key] = genotypes
end
############################################################################
# Read model into environment
############################################################################
M_indices = []
for (i, line) in enumerate(file_data)
    if startswith(line, "#model")
        push!(M_indices, i)
    end
end
model_string = []

for i in eachindex(M_indices)
    current_index = M_indices[i]
    model_sep =  strip(split(file_data[current_index], ':')[2])
    push!(model_string, model_sep)
end

models = "";

# Iterate through the string array, adding each string and a newline character to the merged string.
for str in model_string
    models *= "$str\n"
end

buildmodel_start = findfirst(line -> startswith(line, "#buildmodel"), file_data)
buildmodel_end = findlast(line -> startswith(line, "#buildmodel"), file_data)
file_data_buildmodel = if isnothing(buildmodel_start) || isnothing(buildmodel_end)
    []
else
    file_data[buildmodel_start:buildmodel_end]
end

buildmodel_R_line = findfirst(line -> startswith(line, "#buildmodel_R"), file_data_buildmodel)
buildmodel_df_line = findfirst(line -> startswith(line, "#buildmodel_df"), file_data_buildmodel)
buildmodel_num_hidden_nodes_line = findfirst(line -> startswith(line, "#buildmodel_num_hidden_nodes"), file_data_buildmodel)
buildmodel_nonlinear_function_line = findfirst(line -> startswith(line, "#buildmodel_nonlinear_function"), file_data_buildmodel)
buildmodel_latent_traits_line = findfirst(line -> startswith(line, "#buildmodel_latent_traits"), file_data_buildmodel)
buildmodel_user_σ2_yobs_line = findfirst(line -> startswith(line, "#buildmodel_user_σ2_yobs"), file_data_buildmodel)
buildmodel_user_σ2_weightsNN_line = findfirst(line -> startswith(line, "#buildmodel_user_σ2_weightsNN"), file_data_buildmodel)
buildmodel_censored_trait_line = findfirst(line -> startswith(line, "#buildmodel_censored_trait"), file_data_buildmodel)
buildmodel_categorical_trait_line = findfirst(line -> startswith(line, "#buildmodel_categorical_trait"), file_data_buildmodel)



if buildmodel_R_line !== nothing
    buildmodel_R_value = strip(split(file_data_buildmodel[buildmodel_R_line], ':')[2])
else
    buildmodel_R_value = false
end

if buildmodel_df_line !== nothing
    buildmodel_df_value = strip(split(file_data_buildmodel[buildmodel_df_line], ':')[2])
else
    buildmodel_df_value = 4.0
end

if buildmodel_num_hidden_nodes_line !== nothing
    buildmodel_num_hidden_nodes_value = strip(split(file_data_buildmodel[buildmodel_num_hidden_nodes_line], ':')[2])
else
    buildmodel_num_hidden_nodes_value = false
end

if buildmodel_nonlinear_function_line !== nothing
    buildmodel_nonlinear_function_value = strip(split(file_data_buildmodel[buildmodel_nonlinear_function_line], ':')[2])
else
    buildmodel_nonlinear_function_value = false
end

if buildmodel_latent_traits_line !== nothing
    buildmodel_latent_traits_value = strip(split(file_data_buildmodel[buildmodel_latent_traits_line], ':')[2])
else
    buildmodel_latent_traits_value = false
end

if buildmodel_user_σ2_yobs_line !== nothing
    buildmodel_user_σ2_yobs_value = strip(split(file_data_buildmodel[buildmodel_user_σ2_yobs_line], ':')[2])
else
    buildmodel_user_σ2_yobs_value = false
end

if buildmodel_user_σ2_weightsNN_line !== nothing
    buildmodel_user_σ2_weightsNN_value = strip(split(file_data_buildmodel[buildmodel_user_σ2_weightsNN_line], ':')[2])
else
    buildmodel_user_σ2_weightsNN_value = false
end

if buildmodel_censored_trait_line !== nothing
    buildmodel_censored_trait_value = strip(split(file_data_buildmodel[buildmodel_censored_trait_line], ':')[2])
else
    buildmodel_censored_trait_value = false
end

if buildmodel_categorical_trait_line !== nothing
    buildmodel_categorical_trait_value = strip(split(file_data_buildmodel[buildmodel_categorical_trait_line], ':')[2])
else
    buildmodel_categorical_trait_value = false
end



model = build_model(models,
                    R = buildmodel_R_value,
                    df = buildmodel_df_value, 
                    g_types = g_type, 
                    genoterms = genoterm, 
                    num_hidden_nodes = buildmodel_num_hidden_nodes_value, 
                    nonlinear_function = buildmodel_nonlinear_function_value, 
                    latent_traits=buildmodel_latent_traits_value, #nonlinear_function(x1,x2) = x1+x2
                    user_σ2_yobs = buildmodel_user_σ2_weightsNN_value, 
                    user_σ2_weightsNN = buildmodel_censored_trait_value,
                    censored_trait = buildmodel_censored_trait_value, 
                    categorical_trait = buildmodel_categorical_trait_value);
############################################################################
# set random and covariates 
############################################################################

# Set Factors or Covariates
C_indices = []
for (i, line) in enumerate(file_data)
    if startswith(line, "#covariates")
        push!(C_indices, i)
    end
end

covariates_str = []
for i in eachindex(C_indices)
    current_index = C_indices[i]
    covariates_sep =  strip(split(file_data[current_index], ':')[2])
    push!(covariates_str, covariates_sep)
end

set_covariate(model,covariates_str...);

# Set Random or Fixed Effects without pedigree
randomorfixed_start = findfirst(line -> startswith(line, "#randomorfixed"), file_data)
randomorfixed_end = findlast(line -> startswith(line, "#randomorfixed"), file_data)
file_data_randomorfixed = if isnothing(randomorfixed_start) || isnothing(randomorfixed_end)
    []
else
    file_data[randomorfixed_start:randomorfixed_end]
end

randomorfixed_randomStr_line = findfirst(line -> startswith(line, "#randomorfixed_randomStr"), file_data_randomorfixed)
randomorfixed_G_line = findfirst(line -> startswith(line, "#randomorfixed_G"), file_data_randomorfixed)
randomorfixed_Vinv_line = findfirst(line -> startswith(line, "#randomorfixed_Vinv"), file_data_randomorfixed)
randomorfixed_names_line = findfirst(line -> startswith(line, "#randomorfixed_names"), file_data_randomorfixed)
randomorfixed_df_line = findfirst(line -> startswith(line, "#randomorfixed_df"), file_data_randomorfixed)

if randomorfixed_G_line !== nothing
    randomorfixed_G_value = strip(split(file_data_randomorfixed[randomorfixed_G_line], ':')[2])
else
    randomorfixed_G_value = false
end

if randomorfixed_Vinv_line !== nothing
    randomorfixed_Vinv_value = strip(split(file_data_randomorfixed[randomorfixed_Vinv_line], ':')[2])
else
    randomorfixed_Vinv_value = 0
end

if randomorfixed_names_line !== nothing
    randomorfixed_names_value = strip(split(file_data_randomorfixed[randomorfixed_names_line], ':')[2])
else
    randomorfixed_names_value = []
end

if randomorfixed_df_line !== nothing
    randomorfixed_df_value = strip(split(file_data_randomorfixed[randomorfixed_df_line], ':')[2])
else
    randomorfixed_df_value = 4.0
end

if randomorfixed_randomStr_line !== nothing
    randomorfixed_randomStr_value = strip(split(file_data_randomorfixed[randomorfixed_randomStr_line], ':')[2])
    set_random(model,randomorfixed_randomStr_value, randomorfixed_G_value , Vinv= randomorfixed_Vinv_value,names=randomorfixed_names_value,df=randomorfixed_df_value);
else
    println("NO random effect set in this model")
end

# Set Random or Fixed Effects with pedigree
randomorwithpedigree_start = findfirst(line -> startswith(line, "#randomorwithpedigree"), file_data)
randomorwithpedigree_end = findlast(line -> startswith(line, "#randomorwithpedigree"), file_data)
file_data_randomorwithpedigree = if isnothing(randomorwithpedigree_start) || isnothing(randomorwithpedigree_end)
    []
else
    file_data[randomorwithpedigree_start:randomorwithpedigree_end]
end

randomorwithpedigree_randomStr_line = findfirst(line -> startswith(line, "#randomorwithpedigree_randomStr"), file_data_randomorwithpedigree)
randomorwithpedigree_G_line = findfirst(line -> startswith(line, "#randomorwithpedigree_G"), file_data_randomorwithpedigree)
randomorwithpedigree_df_line = findfirst(line -> startswith(line, "#randomorwithpedigree_df"), file_data_randomorwithpedigree)

if randomorwithpedigree_G_line !== nothing
    randomorwithpedigree_G_value = strip(split(file_data_randomorwithpedigree[randomorwithpedigree_G_line], ':')[2])
else
    randomorwithpedigree_G_value = false
end

if randomorwithpedigree_df_line !== nothing
    randomorwithpedigree_df_value = strip(split(file_data_randomorwithpedigree[randomorwithpedigree_df_line], ':')[2])
else
    randomorwithpedigree_df_value = 4.0
end

if randomorwithpedigree_randomStr_line !== nothing
    randomorwithpedigree_randomStr_value = strip(split(file_data_randomorwithpedigree[randomorwithpedigree_randomStr_line], ':')[2])
    set_random(model, randomorwithpedigree_randomStr_value,pedigree,randomorwithpedigree_G_value,df=randomorwithpedigree_df_value);
else
    println("NO random effect with pedigree set in this model")
end

############################################################################
# Run MCMC analysis
############################################################################

runparameter_start = findfirst(line -> startswith(line, "#runparameter"), file_data)
runparameter_end = findlast(line -> startswith(line, "##runparameter"), file_data)
file_data_runparameter = if isnothing(runparameter_start) || isnothing(runparameter_end)
    []
else
    file_data[runparameter_start:runparameter_end]
end

runparameter_heterogeneous_residuals_line = findfirst(line -> startswith(line, "#runparameter_heterogeneous_residuals"), file_data_runparameter)
runparameter_chain_length_line = findfirst(line -> startswith(line, "#runparameter_chain_length"), file_data_runparameter)
runparameter_starting_value_function_line = findfirst(line -> startswith(line, "#runparameter_starting_value_function"), file_data_runparameter)
runparameter_burnin_line = findfirst(line -> startswith(line, "#runparameter_burnin"), file_data_runparameter)
runparameter_output_samples_frequency_line = findfirst(line -> startswith(line, "#runparameter_output_samples_frequency"), file_data_runparameter)
runparameter_update_priors_frequency_line = findfirst(line -> startswith(line, "#runparameter_update_priors_frequency"), file_data_runparameter)
runparameter_estimate_variance_line = findfirst(line -> startswith(line, "#runparameter_estimate_variance"), file_data_runparameter)
runparameter_single_step_analysis_line = findfirst(line -> startswith(line, "#runparameter_single_step_analysis"), file_data_runparameter)
runparameter_pedigree_line = findfirst(line -> startswith(line, "#runparameter_pedigree"), file_data_runparameter)
runparameter_fitting_J_vector_line = findfirst(line -> startswith(line, "#runparameter_fitting_J_vector"), file_data_runparameter)
runparameter_causal_structure_line = findfirst(line -> startswith(line, "#runparameter_causal_structure"), file_data_runparameter)
runparameter_mega_trait_line = findfirst(line -> startswith(line, "#runparameter_mega_trait"), file_data_runparameter)
runparameter_missing_phenotypes_line = findfirst(line -> startswith(line, "#runparameter_missing_phenotypes"), file_data_runparameter)
runparameter_constraint_line = findfirst(line -> startswith(line, "#runparameter_constraint"), file_data_runparameter)
runparameter_RRM_line = findfirst(line -> startswith(line, "#runparameter_RRM"), file_data_runparameter)
runparameter_outputEBV_line = findfirst(line -> startswith(line, "#runparameter_outputEBV"), file_data_runparameter)
runparameter_output_heritability_line = findfirst(line -> startswith(line, "#runparameter_output_heritability"), file_data_runparameter)
runparameter_prediction_equation_line = findfirst(line -> startswith(line, "#runparameter_prediction_equation"), file_data_runparameter)
runparameter_seed_line = findfirst(line -> startswith(line, "#runparameter_seed"), file_data_runparameter)
runparameter_printout_model_info_line = findfirst(line -> startswith(line, "#runparameter_printout_model"), file_data_runparameter)
runparameter_printout_frequency_line = findfirst(line -> startswith(line, "#runparameter_printout_frequency"), file_data_runparameter)
runparameter_big_memory_line = findfirst(line -> startswith(line, "#runparameter_big_memory"), file_data_runparameter)
runparameter_double_precision_line = findfirst(line -> startswith(line, "#runparameter_double_precision"), file_data_runparameter)
runparameter_output_folder_line = findfirst(line -> startswith(line, "#runparameter_output_folder"), file_data_runparameter)
runparameter_output_samples_for_all_parameters_line = findfirst(line -> startswith(line, "#runparameter_output_samples_for_all_parameters"), file_data_runparameter)
runparameter_methods_line = findfirst(line -> startswith(line, "#runparameter_methods"), file_data_runparameter)
runparameter_Pi_info_line = findfirst(line -> startswith(line, "#runparameter_Pi_info"), file_data_runparameter)
runparameter_estimatePi_line = findfirst(line -> startswith(line, "#runparameter_estimatePi"), file_data_runparameter)
runparameter_estimateScale_line = findfirst(line -> startswith(line, "#runparameter_estimateScale"), file_data_runparameter)
runparameter_categorical_trait_line = findfirst(line -> startswith(line, "#runparameter_categorical_trait"), file_data_runparameter)
runparameter_censored_trait_line = findfirst(line -> startswith(line, "#runparameter_censored_trait"), file_data_runparameter)









if runparameter_heterogeneous_residuals_line !== nothing
    runparameter_heterogeneous_residuals_value = strip(split(file_data_runparameter[runparameter_heterogeneous_residuals_line], ':')[2])
else
    runparameter_heterogeneous_residuals_value = false
end

if runparameter_chain_length_line !== nothing
    runparameter_chain_length_value = strip(split(file_data_runparameter[runparameter_chain_length_line], ':')[2])
else
    runparameter_chain_length_value = 100
end

if runparameter_starting_value_function_line !== nothing
    runparameter_starting_value_function_value = strip(split(file_data_runparameter[runparameter_starting_value_function_line], ':')[2])
else
    runparameter_starting_value_function_value = false
end


if runparameter_burnin_line !== nothing
    runparameter_burnin_value = strip(split(file_data_runparameter[runparameter_burnin_line], ':')[2])
else
    runparameter_burnin_value = 0
end

if runparameter_output_samples_frequency_line !== nothing
    runparameter_output_samples_frequency_value = strip(split(file_data_runparameter[runparameter_output_samples_frequency_line], ':')[2])
else
    runparameter_output_samples_frequency_value = runparameter_chain_length_value>1000 ? div(chain_length,1000) : 1
end

if runparameter_update_priors_frequency_line !== nothing
    runparameter_update_priors_frequency_value = strip(split(file_data_runparameter[runparameter_update_priors_frequency_line], ':')[2])
else
    runparameter_update_priors_frequency_value = 0
end

if runparameter_estimate_variance_line !== nothing
    runparameter_estimate_variance_value = strip(split(file_data_runparameter[runparameter_estimate_variance_line], ':')[2])
else
    runparameter_estimate_variance_value = true
end

if runparameter_single_step_analysis_line !== nothing
    runparameter_single_step_analysis_value = strip(split(file_data_runparameter[runparameter_single_step_analysis_line], ':')[2])
else
    runparameter_single_step_analysis_value = false
end

if runparameter_pedigree_line !== nothing
    runparameter_pedigree_value = strip(split(file_data_runparameter[runparameter_pedigree_line], ':')[2])
else
    runparameter_pedigree_value = false
end

if runparameter_fitting_J_vector_line !== nothing
    runparameter_fitting_J_vector_value = strip(split(file_data_runparameter[runparameter_fitting_J_vector_line], ':')[2])
else
    runparameter_fitting_J_vector_value = true
end

if runparameter_causal_structure_line !== nothing
    runparameter_causal_structure_value = strip(split(file_data_runparameter[runparameter_causal_structure_line], ':')[2])
else
    runparameter_causal_structure_value = false
end

if runparameter_mega_trait_line !== nothing
    runparameter_mega_trait_value = strip(split(file_data_runparameter[runparameter_mega_trait_line], ':')[2])
else
    runparameter_mega_trait_value = buildmodel_nonlinear_function_value == false ? false : true
end

if runparameter_missing_phenotypes_line !== nothing
    runparameter_missing_phenotypes_value = strip(split(file_data_runparameter[runparameter_missing_phenotypes_line], ':')[2])
else
    runparameter_missing_phenotypes_value = buildmodel_nonlinear_function_value == false ? false : true
end

if runparameter_constraint_line !== nothing
    runparameter_constraint_value = strip(split(file_data_runparameter[runparameter_constraint_line], ':')[2])
else
    runparameter_constraint_value = false
end

if runparameter_RRM_line !== nothing
    runparameter_RRM_value = strip(split(file_data_runparameter[runparameter_RRM_line], ':')[2])
else
    runparameter_RRM_value = false
end

if runparameter_outputEBV_line !== nothing
    runparameter_outputEBV_value = strip(split(file_data_runparameter[runparameter_outputEBV_line], ':')[2])
else
    runparameter_outputEBV_value = true
end

if runparameter_output_heritability_line !== nothing
    runparameter_output_heritability_value = strip(split(file_data_runparameter[runparameter_output_heritability_line], ':')[2])
else
    runparameter_output_heritability_value = true
end

if runparameter_prediction_equation_line !== nothing
    runparameter_prediction_equation_value = strip(split(file_data_runparameter[runparameter_prediction_equation_line], ':')[2])
else
    runparameter_prediction_equation_value = false
end

if runparameter_seed_line !== nothing
    runparameter_seed_value = strip(split(file_data_runparameter[runparameter_prediction_equation_line], ':')[2])
else
    runparameter_seed_value = false
end

if runparameter_printout_model_info_line !== nothing
    runparameter_printout_model_info_value = strip(split(file_data_runparameter[runparameter_printout_model_info_line], ':')[2])
else
    runparameter_printout_model_info_value = true
end

if runparameter_printout_frequency_line !== nothing
    runparameter_printout_frequency_value = strip(split(file_data_runparameter[runparameter_printout_frequency_line], ':')[2])
else
    runparameter_printout_frequency_value = runparameter_chain_length_value + 1
end

if runparameter_big_memory_line !== nothing
    runparameter_big_memory_value = strip(split(file_data_runparameter[runparameter_big_memory_line], ':')[2])
else
    runparameter_big_memory_value = false
end

if runparameter_double_precision_line !== nothing
    runparameter_double_precision_value = strip(split(file_data_runparameter[runparameter_double_precision_line], ':')[2])
else
    runparameter_double_precision_value = false
end

if runparameter_output_folder_line !== nothing
    runparameter_output_folder_value = strip(split(file_data_runparameter[runparameter_output_folder_line], ':')[2])
else
    runparameter_output_folder_value = "results"
end

if runparameter_output_samples_for_all_parameters_line !== nothing
    runparameter_output_samples_for_all_parameters_value = strip(split(file_data_runparameter[runparameter_output_samples_for_all_parameters_line], ':')[2])
else
    runparameter_output_samples_for_all_parameters_value = false
end

if runparameter_methods_line !== nothing
    runparameter_methods_value = strip(split(file_data_runparameter[runparameter_methods_line], ':')[2])
else
    runparameter_methods_value = "conventional (no markers)"
end

if runparameter_Pi_info_line !== nothing
    runparameter_Pi_info_value = strip(split(file_data_runparameter[runparameter_Pi_info_line], ':')[2])
else
    runparameter_Pi_info_value = 0.0
end

if runparameter_estimatePi_line !== nothing
    runparameter_estimatePi_value = strip(split(file_data_runparameter[runparameter_estimatePi_line], ':')[2])
else
    runparameter_estimatePi_value = false
end

if runparameter_estimateScale_line !== nothing
    runparameter_estimateScale_value = strip(split(file_data_runparameter[runparameter_estimateScale_line], ':')[2])
else
    runparameter_estimateScale_value = false
end

if runparameter_categorical_trait_line !== nothing
    runparameter_categorical_trait_value = strip(split(file_data_runparameter[runparameter_categorical_trait_line], ':')[2])
else
    runparameter_categorical_trait_value = false
end

if runparameter_censored_trait_line !== nothing
    runparameter_censored_trait_value = strip(split(file_data_runparameter[runparameter_censored_trait_line], ':')[2])
else
    runparameter_censored_trait_value = false
end

out=runMCMC(model,phenotypes,
heterogeneous_residuals           = runparameter_heterogeneous_residuals_value,
#MCMC
chain_length                      = runparameter_chain_length_value,##important
starting_value                    = runparameter_starting_value_function_value,
burnin                            = runparameter_burnin_value,##important
output_samples_frequency          = runparameter_output_samples_frequency_value,##important
update_priors_frequency           = runparameter_update_priors_frequency_value,
#Methods
estimate_variance               = runparameter_estimate_variance_value,  ##important
single_step_analysis            = runparameter_single_step_analysis_value, #parameters for single-step analysis  ##important
pedigree                        = runparameter_pedigree_value, #parameters for single-step analysis  ##important

fitting_J_vector                = runparameter_fitting_J_vector_value,  #parameters for single-step analysis
causal_structure                = runparameter_causal_structure_value,
mega_trait                      = runparameter_mega_trait_value, #NNBayes -> default mega_trait=true
missing_phenotypes              = runparameter_missing_phenotypes_value, #NN-MM -> missing hidden nodes will be sampled
constraint                      = runparameter_constraint_value,
RRM                             = runparameter_RRM_value, #  RRM -> false or a matrix (Phi: orthogonalized time covariates)
#Genomic Prediction
outputEBV                       = runparameter_outputEBV_value,
output_heritability             = runparameter_output_heritability_value,
prediction_equation             = runparameter_prediction_equation_value,
#MISC
seed                            = runparameter_seed_value,
printout_model_info             = runparameter_printout_model_info_value,
printout_frequency              = runparameter_printout_frequency_value,
big_memory                      = runparameter_big_memory_value,
double_precision                = runparameter_double_precision_value,
#MCMC samples (defaut to marker effects and hyperparametes (variance componets))
output_folder                     = runparameter_output_folder_value, ##important
output_samples_for_all_parameters = runparameter_output_samples_for_all_parameters_value,
#for deprecated JWAS
methods                         = runparameter_methods_value,
Pi                              = runparameter_Pi_info_value,
estimatePi                      = runparameter_estimatePi_value,
estimateScale                   = runparameter_estimateScale_value,
categorical_trait               = runparameter_categorical_trait_value,  #this has been moved to build_model()
censored_trait                  = runparameter_censored_trait_value) ;

############################################################################
# Check Accuracy
############################################################################





end

function julia_main()
    try
        readapp()
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end
end
