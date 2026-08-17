#=
    Script to fix categories missing and stuff
=#

using DataFrames, CSV

#=
    Functions
=#

#=
    Start the script
=#

# Paths and directories
DIR_MAIN = pwd();
DIR_DATA = joinpath(DIR_MAIN,"analyzed_data");
FILE_DAT = "dat.csv";
DIR_SAVE = joinpath(DIR_MAIN,"figures");

# Read the directory 
files=readdir(DIR_DATA);

# Get only those of the structure factor
files=filter(s -> occursin("structure_factor_", s), files);

# Read the files
df_files=[CSV.read(joinpath(DIR_DATA,file), DataFrame) for file in files];

# The idea is that all df have the exact sames keys

keys_per_simulation = names.(df_files);

# Get the unique keys
keys = unique(keys_per_simulation);

if length(keys) == 1
    println("Same amount of keys")
end

# Get all the keys
all_keys = union(keys...);

# Amount of keys. For first test of similarity
n_keys = length(all_keys);

# Find the simulation that does not have the same keys
diff_keys_per_simulation = map(k_p_s->mapreduce(s->k_p_s.== s,+,all_keys),keys_per_simulation);

# Known wich are equal and what are different
mask_diff = sum.(diff_keys_per_simulation) .!= n_keys;

# Get the files that have different keys
files_diff = files[mask_diff];
