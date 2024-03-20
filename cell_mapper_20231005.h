#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
// #include <float.h>

struct BCarray{
    int zone;
    int type;
    double value;
};

int uniquenode(char *,char *,struct BCarray *,int,int,int ,double ,double ,double ,double);
void circle_path_generator(int** circle_path,int);
int ply_parser(double **nodes,int **cells,FILE *ifp,int num_nodes,char *line,int max_line_size);
int cmp(const void  *a, const void *b);
int edge_comparison(int** edge_to_cell,int** edges,int** boundary_edges,int** cells,int** node_list,int npc,int* num_edges,int* num_boundary_edges,int skip,int edge_count);
void freeMatrix_i2(int **matrix ,int row);
void freeMatrix_d2(double **matrix ,int row);
void freeMatrix_d3(double ***matrix ,int row,int col);
double *** allocateMatrix_d3(int row,int col,int elems);
int ** allocateMatrix_i2(int row,int elems);
double ** allocateMatrix_d2(int row,int elems);
double ** allocateMatrix_d2_zeros(int row,int elems);
void node_weights(double** weight_node,double* weight_node_sum,double** nodes,int** cells,double** cell_center,int num_cells,int npc);
void cell_quantities(int** circle_path,double*** vertices,int** cells,double** nodes,double*** t,double*** t̂,double*** face_mid,double** n_cell,double** n̂_cell,double*** n̂_face,double* volume_cell,double** d_cell_face,double** cell_center,int num_cells,int npc,int vec_size);
double dot_3(double* a,double* b);
double norm_3(double* a);
void cross_product_3(double* a, double α_a, double* b, double α_b, double* c);
void vector_add_3(double* a, double α_a, double* b, double α_b, double* c);
void vector_mult_3(double* a, double α_a, double* c);
void edge_quantities(int** edges,double** nodes,double** d_cell_face,double* area_face,double* area_bface,double** cell_center,double*** n̂_face,double*** face_mid,double** weight_face,double* weight_bface,double*** lf,double** lbf,double* δ_face,double* δ_bface,int num_edges,int vec_size);
void IC(int** cells,double* phi,double ICvalue,int num_cells);
void BC(int** boundary_edges,double* phi_bface,double* phi_node,struct BCarray* BCzone,int num_boundary_edges,int num_boundary_conditions);
void grid_properties(int** cells, int** edges,double* aₚ,double** aₙ,double* sourceᵤ,double* sourceₚ,double* sourceₛ,double** t̂_dot_lf,double*** t̂,double*** lf,double* area_face,double* δ_face,double Gamma,double Sourceᵤ,double Sourceₚ,double* volume_cell, struct BCarray* BCzone,int num_cells,int npc,int num_boundary_conditions);
void gauss_seidel(int** circle_path,int iterations,int** cells,double* phi_node,double* phi,double** weight_node,double* δ_face,double Gamma,int** boundary_edges,struct BCarray* BCzone,double* aₚ,double** aₙ,double* sourceᵤ,double* sourceₚ,double* sourceₛ,double** t̂_dot_lf,int npc,int num_cells,int num_boundary_edges,int num_boundary_conditions);
void bubbleSort(int arr[], int n);
void swap(int* arr, int i, int j);

