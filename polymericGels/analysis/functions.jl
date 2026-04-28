"""
    File with all the functions needed for multiple.jl script
"""
#######
#   STRUCTURE FACTOR RELATED FUNCTIONS
#######

function getDump(dir,file_name)
"""
    Get the data from a single dump file that stores one timeste information
"""
    data = split.(readlines(joinpath(dir,file_name))," ")[9:end];
    HEADERS=data[1][3:end];
    INFO=parse.(Float64,reduce(hcat,data[2:end]))';

    return DataFrame(INFO,HEADERS)
end

function dotSpherical(theta,phi,r)
"""
    Compute the dot product betwen a position and a unit vector r in psherical coordinates.
"""
    q_x=cos(theta)*sin(phi);
    q_y=sin(theta)*sin(phi);
    q_z=cos(phi);
    return q_x*r[1]+q_y*r[2]+q_z*r[3]
end

function densityRhoQ(q_mag,dot_qr)
"""
    Compute the squared of the absolute value of the density at the reciprocal space.
    |rho(r)|^2 = A(vec{q}cdotvec{r})^2 + B(vec{q}cdotvec{r}^2)
    A = sumcos(); B = sumsin()
"""
    return sum(cos.(q_mag*dot_qr))^2 + sum(sin.(q_mag*dot_qr))^2
end

function computeDensity(theta,phi,lambda_o,lambda_f,N_lambda,r)
"""
    Function that computes the static structure factor de different wave vectors.
    Returns a vector with the following interpretation of the values:
    [row] -> [magnitude]}
"""

    # Calculo del producto punto
    dot_qr=[dotSpherical(theta[s],phi[s],r) for s in 1:length(theta)];

    # Evaluación de la densidad y promedio
    # [renglon x columna] -> [ mag x direccion ]
    q_min=2*pi/lambda_f;
    q_max=2*pi/lambda_o;
    q_dom=range(q_min,q_max,length=N_lambda);

    rho_q=[densityRhoQ(l,d) for l in q_dom, d in dot_qr];

    return [q_dom,reduce(vcat,mean(rho_q,dims=2))]

end

function getPosition(dump)
"""
    Get the position of the central particles of a given dump
"""

    # Filtrar
    mask=(dump.type .==1) .| (dump.type .== 2.0);
    dump_filtered=dump[mask,:];

   return [dump_filtered.x,dump_filtered.y,dump_filtered.z]

end

function structureFactor(theta,phi,lambda_o,lambda_f,N_lambda,r_exp)
"""
    Compute the static structure factor
"""
    data=[computeDensity(theta,phi,lambda_o,lambda_f,N_lambda,r) for r in r_exp];
    data=reduce(hcat,data);
    q_domain=collect(first(unique(data[1,:])));
    Sq=reduce(vcat,mean(reduce(hcat,data[2,:]),dims=2));

    return [q_domain,Sq]
end

function getTimeEvolSq(N_qu,dump_paths,time_instant)
"""
    Compute the time evolution of the structure factor
"""

    # Vector unitario del vector de onda
    theta=2*pi*rand(N_qu);
    phi=pi*rand(N_qu); 

    # Obtenemos los dumps de los N experimentos para un instante de tiempo
    dumps=[getDump(path,time_instant) for path in dump_paths];

    r_exp=[getPosition(df) for df in dumps];

    Sq=structureFactor(theta,phi,lambda_o,lambda_f,N_lambda,r_exp)
   
    return Sq

end

#######
#   OTHER  FUNCTIONS
#######

function meanFixystem(dirs)
"""
    Function that returns a dataframe with the mean of N experiments of the observables stored in a fix file.
"""

    # Obtener la información de lo fix files
    data_fix=map(s->dataSystem=extractFixScalar(s,"system_assembly.fixf"),dirs);

    # Variables auxiliares
    n_row = nrow(data_fix[1]);
    n_col = ncol(data_fix[1]);
    cols = names(data_fix[1]);

    # Apilar matrices en un arreglo 3D
    arr = Array{Float64}(undef, n_row, n_col, length(data_fix))
    for (i, df) in enumerate(data_fix)
        arr[:, :, i] = Matrix(df)   # asume columnas numéricas
    end

    # Promedio a lo largo de la tercera dimensión
    promedio = dropdims(mean(arr, dims=3), dims=3);

    return DataFrame(promedio,cols);
end


function extractFixScalar(path_system,file_name)
"""
    Function that extracts the information of fix files that stores global scalar values
"""
    aux=split.(readlines(joinpath(path_system,file_name))," ");
    header=aux[2][2:end];
    info=reduce(hcat,map(s->parse.(Float64,s),aux[3:end]));

    return DataFrame(info',header)
end

function getDat(path)
"""
    Creates a dataframe from the dat file of the experiment.
"""
    file_path=joinpath(path,"dataAssembly.dat");
    aux=split.(readlines(file_path),",");

    df_aux=DataFrame();
    for (col, val) in zip(aux[1], aux[2])
        # Convertir a Float64 si es posible, sino mantener como String
        parsed_val = tryparse(Float64, val)
        if parsed_val !== nothing
            df_aux[!, col] = [parsed_val]
        else
            df_aux[!, col] = [val]
        end
    end

    df_aux[!, :dir] = [last(split(path,"/"))];

    return df_aux
end


