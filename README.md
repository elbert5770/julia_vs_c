# julia_vs_c
File manipulations, allocations, Gauss-Seidel compared between Julia and C.

A single threaded task on Intel(R) Xeon(R) Gold 6230 CPU @ 2.10GHz

**C**
gcc -O3 -o cell_mapper cell_mapper_20240319.c -lm

the full code runs in 8.53 seconds.

gcc -ggdb -o cell_mapper cell_mapper_20231006.c -lm

the full code runs in 19.0 seconds.


**Julia**

julia> include("cell_locations_20240319_msh_ply_3D_views.jl")

@time records: 8.47 seconds

julia> include("cell_locations_20231013_msh_ply_3D_noinb.jl")

13.3 seconds

The conclusion is that with similar levels of optimization

