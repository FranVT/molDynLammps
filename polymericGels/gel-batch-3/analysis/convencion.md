# Para marcar las convenciones de sintaxis

## Definicion de variables
- Variables van en usan minúsculas (nombre)
- Variables se usa guión bajo para separar (apellido)
- Constantes globales y tipos de datos se emplea CamelCase
- Constantes numéricas en mayúsculas y guiones bajos

```julia
aux_dir=joinpath(pwd,df.dir[2]);
Npart=df.Npart[2];
G=9.81;
```


## Comentario de funciones
```julia
"""
    calcular_area(radio)

Devuelve el área de un círculo dado su `radio`.

# Argumentos
- `radio::Real`: El radio del círculo, que debe ser un número positivo.

# Ejemplos
```julia-repl
julia> calcular_area(2.0)
12.566370614359172
"""
function calcular_area(radio)
    return π * radio^2
end
```
