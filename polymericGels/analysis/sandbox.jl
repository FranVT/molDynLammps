# Get directories 
MAIN_DIR=pwd();
DAT_PATH=joinpath(MAIN_DIR,"datFiles","experiments_dat.csv");
SAVE_DIR=joinpath(MAIN_DIR,"analyzedData");

dat_files=CSV.read(DAT_PATH,DataFrame);

# Selection of the system by parameters
phi=0.05;
Temp=0.05;
N_part=5000.0;
CL_con=0.05;

# Se filtra el dataframe 
dat_DF = subset(dat_files,
    :phi => ByRow(==(phi)),
    :Temperature => ByRow(==(Temp)),
    :Npart => ByRow(==(N_part)),
    :"CL-Con" => ByRow(==(CL_con))
)

# Path to the dumps
dump_paths=joinpath.(dat_DF.PARENT_DIR,dat_DF.dir,"traj");

# Parametros para obtener el factor de estructura
N_qu=2^7; # EXPONENTE DEBE SER IMPAR Cantidad de direcciones
lambda_o=0.5; # Limites del rango a explorar (Monomero)
lambda_f=2*dat_DF.L[1]; # Limites del rango a explorar (Tamaño de la caja)
N_lambda=2^9; # Cantidad de magnitudes
N_instants=2;

# Seleccion de time instants
aux_timeStep=Int.((0:dat_DF."save-dump"[1]:(dat_DF."N_heat"[1] + dat_DF."N_isot"[1])));
ind=round.(Int, LinRange(1, length(aux_timeStep), N_instants));
aux_id=aux_timeStep[ind];

time_instants=[replace("traj_assembly.*.dumpf", "*" => string(it)) for it in aux_id];



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


N_qu=2^7; # EXPONENTE DEBE SER IMPAR Cantidad de direcciones
lambda_o=0.5; # Limites del rango a explorar (Monomero)
lambda_f=2*dat_DF.L[1]; # Limites del rango a explorar (Tamaño de la caja)
N_lambda=2^9; # Cantidad de magnitudes
N_instants=2;

time_instant=time_instants[1];

    # Vector unitario del vector de onda
    N_phi=Int64(sqrt(div(N_qu,2)));
    N_theta=Int64(2*N_phi);

    theta=2*pi*rand(N_theta);
    phi=pi*rand(N_phi); 

    # Obtenemos los dumps de los N experimentos para un instante de tiempo
    dumps=[getDump(path,time_instant) for path in dump_paths];

    r_exp=[getPosition(df) for df in dumps];

    r=r_exp[1];


    # Calculo del producto punto
    dot_qr=[dotSpherical(th,ph,r) for th in theta, ph in phi];

    # Evaluación de la densidad y promedio
    # [renglon x columna] -> [ mag x direccion ]
    q_min=2*pi/lambda_f;
    q_max=2*pi/lambda_o;
    q_dom=range(q_min,q_max,length=N_lambda);

    rho_q=[densityRhoQ(l,d) for l in q_dom, d in dot_qr];
    rho_q=reduce(vcat,mean(mean(rho_q,dims=3),dims=2))./length(r[1]);


    # Compute the average with the same magnitude, different directions
    output= [q_dom,rho_q];











#=
Sq_t=[getTimeEvolSq(N_qu,dump_paths,time_instant) for time_instant in time_instants];
function getTimeEvolSq(N_qu,dump_paths,time_instant)
"""
    Compute the time evolution of the structure factor
"""

    # Vector unitario del vector de onda
    N_phi=Int64(sqrt(div(N_qu,2)));
    N_theta=Int64(2*N_phi);

    theta=2*pi*rand(N_theta);
    phi=pi*rand(N_phi); 

    # Obtenemos los dumps de los N experimentos para un instante de tiempo
    dumps=[getDump(path,time_instant) for path in dump_paths];

    r_exp=[getPosition(df) for df in dumps];

    Sq=structureFactor(theta,phi,lambda_o,lambda_f,N_lambda,r_exp)
 function structureFactor(theta,phi,lambda_o,lambda_f,N_lambda,r_exp)
"""
    Compute the static structure factor
"""
    data=[computeDensity(theta,phi,lambda_o,lambda_f,N_lambda,r) for r in r_exp];
    function computeDensity(theta,phi,lambda_o,lambda_f,N_lambda,r)
"""
    Function that computes the static structure factor de different wave vectors.
    Returns a vector with the following interpretation of the values:
    [row] -> [magnitude]}
"""

    # Calculo del producto punto
    dot_qr=[dotSpherical(th,ph,r) for th in theta, ph in phi];
function dotSpherical(theta,phi,r)
"""
    Compute the dot product betwen a position and a unit vector r in psherical coordinates.
"""
    q_x=cos(theta)*sin(phi);
    q_y=sin(theta)*sin(phi);
    q_z=cos(phi);
    return q_x*r[1]+q_y*r[2]+q_z*r[3]
end


    # Evaluación de la densidad y promedio
    # [renglon x columna] -> [ mag x direccion ]
    q_min=2*pi/lambda_f;
    q_max=2*pi/lambda_o;
    q_dom=range(q_min,q_max,length=N_lambda);

    rho_q=[densityRhoQ(l,d) for l in q_dom, d in dot_qr];
function densityRhoQ(q_mag,dot_qr)
"""
    Compute the squared of the absolute value of the density at the reciprocal space.
    |rho(r)|^2 = A(vec{q}cdotvec{r})^2 + B(vec{q}cdotvec{r}^2)
    A = sumcos(); B = sumsin()
"""
    return sum(cos.(q_mag*dot_qr))^2 + sum(sin.(q_mag*dot_qr))^2
end


    # Compute the average with the same magnitude, different directions
    return [q_dom,reduce(vcat,mean(rho_q,dims=2))]

end


    data=reduce(hcat,data);
    q_domain=collect(first(unique(data[1,:])));

    # Compute assembly average 
    Sq=reduce(vcat,mean(reduce(hcat,data[2,:]),dims=2))./length(r_exp);

    return [q_domain,Sq]
end

  
    return Sq

end






=#
