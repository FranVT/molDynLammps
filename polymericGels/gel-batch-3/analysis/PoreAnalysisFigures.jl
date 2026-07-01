#=
    Script to create figures of pore distribution
=#

using DataFrames, CSV
using GLMakie, LaTeXStrings 
#=
    Set figures 
=#
INCH = 96;
PT = 4/3;
CM = INCH / 2.54;

set_theme!(
    backgroundcolor = :white,
    fontsize = 16PT,                     # tamaño de letra base
    Axis = (
        titlesize = 20PT,                 # tamaño del título del eje
        xlabelsize = 14PT,                # tamaño etiqueta eje x
        ylabelsize = 14PT,                 # tamaño etiqueta eje y
        xticklabelsize = 12PT,             # tamaño números eje x
        yticklabelsize = 12PT,              # tamaño números eje y
        xgridstyle = :dash,               # estilo de la cuadrícula
        ygridstyle = :dash,
        spinewidth = 1.5PT,
    ),
    Legend = (
        labelsize = 10PT,                   # tamaño texto leyenda
        framewidth = 1.5PT,
    ),
    Colorbar = (
        labelsize = 14PT,
        ticklabelsize = 12PT,
    ),
    Figure = (
        size = (15CM, 12CM)                 # tamaño de la figura (ancho, alto)
    )
)

#=
    Functions 
=#

#=
    Start the script
=#

# Paths and directories
DIR_MAIN = pwd();
DIR_SAVE = joinpath(DIR_MAIN,"analyzed_data");
FILE_DAT = "dat.csv";

# Constants
DT=0.001;

# Read the directory 
files=readdir(DIR_SAVE);

# Get only those of the structure factor
files=filter(s -> occursin("pore_analysis_og_", s), files);

# Read the files
df_files=[CSV.read(joinpath(DIR_SAVE,file), DataFrame) for file in files];


