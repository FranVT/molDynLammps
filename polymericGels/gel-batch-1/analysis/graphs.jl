"""
    Script to make graphs
"""

using CSV, DataFrames

# Get directories 
MAIN_DIR=dirname(pwd());
INFO_DIR=joinpath(pwd(),"data_mod");

DF_DIR=filter(isfile,readdir(INFO_DIR,join=true));

df=map(s->CSV.read(DF_DIR,DataFrame),filter(str -> occursin("dat-", str), DF_DIR));
