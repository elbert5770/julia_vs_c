#include "cell_mapper_20231005.h"


void main(){
   clock_t start_time = clock();
    int num_boundary_conditions = 4;
    struct BCarray *BCzone = malloc(sizeof(struct BCarray) * num_boundary_conditions);
  
    BCzone[0] = (struct BCarray){.zone = 5,.type = 1, .value = 1.0};
    BCzone[1] = (struct BCarray){.zone = 6,.type = 1, .value = 0.0};
    BCzone[2] = (struct BCarray){.zone = 7,.type = 2, .value = 0.0};
    BCzone[3] = (struct BCarray){.zone = 8,.type = 1, .value = 0.0};

    double Sourceᵤ = 0.0;
    double Sourceₚ = 0.0;

    double ICvalue = 0.0;

    double Gamma = 0.1;
    
    int vec_size = 3;
    int npc = 3; // Nodes per cell, 3 or 4
    if (npc != 3 && npc != 4){
      printf("npc (number of cells per node) must be 3 or 4\n");
      return;
    }
    int result;
    
    char cell_number[] = ("864691136812081779");
    // printf("%s\n",cell_number);
    char file_type[] = (".ply");
    char filename[strlen(cell_number)+strlen(file_type)];
    strcpy(filename,cell_number);
    // printf("%s\n",filename);
    strcat(filename,file_type);
    printf("Input file: %s\n",filename);
    result = uniquenode(filename,cell_number,BCzone,num_boundary_conditions,npc,vec_size,ICvalue,Gamma,Sourceᵤ,Sourceₚ);//, &cell_number,  &BC, ICvalue, Gamma, Source_u, Source_p);
    printf("Result code: \n",result);
    
    free(BCzone);
    double elapsed_time = (double)(clock() - start_time) / CLOCKS_PER_SEC;
  printf("Done in %f seconds\n", elapsed_time);
    return;
}
//     # filename = ("C:\\Users\\elber\\Documents\\git\\MicronsBinder\\notebooks\\intro\\864691136812081779_2.ply")
//     # filename = ("C:\\Users\\elber\\Documents\\git\\MicronsBinder\\notebooks\\intro\\synapse_test_2D.msh")
//     # filename = ("C:\\Users\\elber\\Documents\\git\\MicronsBinder\\notebooks\\intro\\fluent_test_filesuser_filesv232file.msh")
//     # filename = ("C:\\Users\\elber\\Documents\\git\\MicronsBinder\\notebooks\\intro\\skewed_diffusion_2D.msh")
//     # filename = ("C:\\Users\\elber\\Documents\\git\\MicronsBinder\\notebooks\\intro\\small_mesh2D_2.msh")
//     # filename = ("C:\\Users\\elber\\Documents\\git\\MicronsBinder\\notebooks\\intro\\rectangular_mesh2D_3.msh")
//     # filename = ("C:\\Users\\elber\\Documents\\git\\MicronsBinder\\notebooks\\intro\\cylinder_thin_tri.msh")
//     # filename = ("C:\\Users\\elber\\Documents\\git\\MicronsBinder\\notebooks\\intro\\synapse_test_2D.msh")
//     filename = ("C:\\Users\\elber\\Documents\\git\\MicronsBinder\\notebooks\\intro\\864691136812081779.ply")
    
//     # filename = ("C:\\Users\\elber\\Documents\\git\\MicronsBinder\\notebooks\\intro\\864691136334080435.ply")

int uniquenode(char *filename,char *cell_number,struct BCarray *BCzone, int num_boundary_conditions,int npc,int vec_size,double ICvalue, double Gamma, double Sourceᵤ, double Sourceₚ){
// printf("Boundary zones\n");
// for (int i = 0; i < 4; i++) {
//   printf("Zone: %d, Type: %d, Value: %f\n", BCzone[i].zone, BCzone[i].type, BCzone[i].value);
// }
// //   printf("%s\n",filename);
// printf("\n",filename);
  int num_nodes = 0;
  int num_cells = 0;

  int num_dim = 3;
  int count = 0,i,j,index;
  FILE *ifp;
  ifp = fopen(filename,"r");
  if (ifp == NULL) {
    // Handle error.
    return 1;
  }

  // Read the file line by line.
  // char line[1024];
  char *line;
  int max_line_size = 1024;
  line = malloc(sizeof(char) * max_line_size);

    // printf("%d\n",sizeof(line));
  while (fgets(line, max_line_size, ifp) != NULL) {
    if (strncmp(line, "element vertex ", strlen("element vertex ")) == 0) {
      num_nodes = strtol(line + 15, NULL, 10);
      printf("Number of nodes: %d\n",num_nodes);
    }

    if (strncmp(line, "element face ", strlen("element face ")) == 0) {
      num_cells = strtol(line + 13, NULL, 0);
      printf("Number of cells: %d\n",num_cells);
    }
    count += 1;
    if (strncmp(line, "end_header", strlen("end_header")) == 0) {
      // printf("%d\n",count);
      break;
    }
    
   
    


  }
  int **circle_path = allocateMatrix_i2(npc,2);
  circle_path_generator(circle_path,npc);
  
  

  double **nodes = allocateMatrix_d2(num_nodes,3);
 

  int cell_length = npc*3+1;
  int **cells = allocateMatrix_i2(num_cells,cell_length);
  

  ply_parser(nodes,cells,ifp,num_nodes,line,max_line_size);
  
  int edge_count = num_cells*npc;
  int **edge_to_cell = allocateMatrix_i2(edge_count,5);


  
  for (i = 0;i < num_cells;i++){
    for (j = 0;j<npc;j++){
        index = i*npc+j;
        // printf("%d %d %d\n",i,j,index);
        if (cells[i][circle_path[j][0]] < cells[i][circle_path[j][1]]){
            
            edge_to_cell[index][0] = cells[i][circle_path[j][0]];
            edge_to_cell[index][1] = cells[i][circle_path[j][1]];
            edge_to_cell[index][4] = 0;
        }else{
            edge_to_cell[index][0] = cells[i][circle_path[j][1]];
            edge_to_cell[index][1] = cells[i][circle_path[j][0]];
            edge_to_cell[index][4] = 1;
        }
        edge_to_cell[index][2] = i;
        edge_to_cell[index][3] = j;           
    }
  }

  qsort(edge_to_cell, edge_count, sizeof(*edge_to_cell), cmp);
 
  int **node_list = allocateMatrix_i2(2,5);

  int skip=1;

  int **edges = allocateMatrix_i2(edge_count,11);
  
  int **boundary_edges = allocateMatrix_i2(edge_count,11);
  
  int num_edges,num_boundary_edges;
  edge_comparison(edge_to_cell,edges,boundary_edges,cells,node_list,npc,&num_edges,&num_boundary_edges,skip,edge_count);
  
  for (i = num_edges;i < edge_count;i++){
    free(edges[i]);
  }
  edges = realloc(edges,sizeof(int *)*num_edges);
  for (i = num_boundary_edges;i < edge_count;i++){
    free(boundary_edges[i]);
  }
  boundary_edges = realloc(boundary_edges,sizeof(int *)*num_boundary_edges);

  double ***vertices = allocateMatrix_d3(num_cells,npc,3);
 

  double ***t = allocateMatrix_d3(num_cells,npc,3);
  double ***t̂ = allocateMatrix_d3(num_cells,npc,3);
  double ***face_mid = allocateMatrix_d3(num_cells,npc,3);
  double ***n̂_face = allocateMatrix_d3(num_cells,npc,3);
  double **n_cell = allocateMatrix_d2(num_cells,3);
  double **n̂_cell = allocateMatrix_d2_zeros(num_cells,3);
  double *volume_cell = malloc(sizeof(double) * num_cells);
  double **d_cell_face = allocateMatrix_d2(num_cells,npc);
  double **cell_center = allocateMatrix_d2_zeros(num_cells,3);
 
  cell_quantities(circle_path,vertices,cells,nodes,t,t̂,face_mid,n_cell,n̂_cell,n̂_face,volume_cell,d_cell_face,cell_center,num_cells,npc,vec_size);

  double *area_face = malloc(sizeof(double) * num_edges);
  double *area_bface = malloc(sizeof(double) * num_boundary_edges);

  double ***lf = allocateMatrix_d3(num_edges,2,3);
  double *δ_face = malloc(sizeof(double) * num_edges);
  double **weight_face = allocateMatrix_d2(num_edges,2);
  double **lbf = allocateMatrix_d2(num_boundary_edges,3);
  double *δ_bface = malloc(sizeof(double) * num_boundary_edges);
  double *weight_bface = malloc(sizeof(double) * num_boundary_edges);


  edge_quantities(edges,nodes,d_cell_face,area_face,area_bface,cell_center,n̂_face,face_mid,weight_face,weight_bface,lf,lbf,δ_face,δ_bface,num_edges,vec_size);
  
  double **weight_node = allocateMatrix_d2(num_cells,npc);
  double *weight_node_sum = calloc(num_nodes,sizeof(double));

  node_weights(weight_node,weight_node_sum,nodes,cells,cell_center,num_cells,npc);

  double *phi = malloc(sizeof(double) * num_cells);
  IC(cells,phi,ICvalue,num_cells);
  double *phi_bface = calloc(num_boundary_edges,sizeof(double));
  double *phi_node = calloc(num_nodes,sizeof(double));


  double **t̂_dot_lf = allocateMatrix_d2(num_cells,npc);

  double *aₚ = calloc(num_cells,sizeof(double));
  double **aₙ = allocateMatrix_d2(num_cells,npc);

  double *sourceᵤ = calloc(num_cells,sizeof(double));
  double *sourceₚ = calloc(num_cells,sizeof(double));
  double *sourceₛ = calloc(num_cells,sizeof(double));

  BC(boundary_edges,phi_bface,phi_node,BCzone,num_boundary_edges,num_boundary_conditions);

  grid_properties(cells, edges,aₚ,aₙ,sourceᵤ,sourceₚ,sourceₛ,t̂_dot_lf,t̂,lf,area_face,δ_face,Gamma,Sourceᵤ,Sourceₚ,volume_cell,BCzone,num_cells,npc,num_boundary_conditions);
 
 int iterations = 1000;
  gauss_seidel(circle_path, iterations, cells, phi_node, phi, weight_node, δ_face, Gamma, boundary_edges, BCzone, aₚ, aₙ, sourceᵤ, sourceₚ, sourceₛ, t̂_dot_lf, npc, num_cells, num_boundary_edges, num_boundary_conditions);


  free(sourceᵤ);
  free(sourceₚ);
  free(sourceₛ);
  freeMatrix_d2(aₙ,num_cells);
  free(aₚ);
  freeMatrix_d2(t̂_dot_lf,num_cells);
  free(phi_node);
  free(phi_bface);
  free(phi);
  freeMatrix_d2(weight_node,num_cells);
  free(weight_node_sum);
  free(weight_bface);
  free(δ_bface);
  freeMatrix_d2(lbf,num_boundary_edges);
  freeMatrix_d2(weight_face ,num_edges);
  free(δ_face);
  freeMatrix_d3(lf,num_edges,2);
  free(area_bface);
  free(area_face);
  freeMatrix_d3(face_mid,num_cells,npc);
  freeMatrix_d3(n̂_face,num_cells,npc);
  freeMatrix_d2(n_cell ,num_cells);
  freeMatrix_d2(n̂_cell ,num_cells);
  freeMatrix_d2(d_cell_face ,num_cells);
  freeMatrix_d2(cell_center ,num_cells);
  free(volume_cell);


 
  free(line);
  freeMatrix_d2(nodes ,num_nodes);

  freeMatrix_i2(cells ,num_cells);


  freeMatrix_i2(edge_to_cell ,edge_count);

      
  freeMatrix_i2(node_list ,2);


    
  freeMatrix_d3(vertices,num_cells,npc);
  freeMatrix_d3(t,num_cells,npc);
  freeMatrix_d3(t̂,num_cells,npc);


  freeMatrix_i2(edges ,num_edges);

  
  freeMatrix_i2(boundary_edges ,num_boundary_edges);
  freeMatrix_i2(circle_path ,npc);
  


  fclose(ifp);

  return 0;
}

void circle_path_generator(int** circle_path,int npc){
  for (int i = 0; i < npc; i++){
      for (int j = 0; j < 2; j++){
        circle_path[i][j] = i+j;
      }
  }

  circle_path[npc-1][1] = 0;
   

}

int ply_parser(double **nodes,int **cells,FILE *ifp,int num_nodes,char *line,int max_line_size){
  char *found;
  int count2 = 0;
  int count3 = 0;
  int count4;
  int i;
  int cellfaces[3];
  int temp;
  // printf("%d\n",sizeof(line));
  // printf("%s",line);
 
  while (fgets(line, max_line_size, ifp) != NULL) {
    
    if (count2 < num_nodes){
      // printf("%d %s",count2,line);
      // nodes[count2-count-num_nodes] = strtod(line);
      found = strtok(line," ");
      nodes[count2][0] = strtod(found,NULL);
      count3 += 1;
      if( found==NULL)
      {
        //   printf("\t'%s'\n",line);
          puts("\tNo separators found");
          return(1);
      }
      while(found && count3 < 3)
      {
        //   printf("\t'%s'\n",found);
          
          found = strtok(NULL," ");
          nodes[count2][count3] = strtod(found,NULL);
          count3 += 1;
      }
      // printf("(%d, %f, %f, %f)\n\n", count2,nodes[count2][0], nodes[count2][1], nodes[count2][2]);
      count2 += 1;
      count3 = 0;
      
    }
    else{
      count4 = count2-num_nodes;

      found = strtok(line," ");
      
      if( found==NULL)
      {
        //   printf("\t'%s'\n",line);
          puts("\tNo separators found");
          return(1);
      }
      while(found && count3 < 3)
      {
              
        found = strtok(NULL," ");
        
        cellfaces[count3]= atoi(found);
        
        count3 += 1;
      }

      bubbleSort(cellfaces,3);
      for (i=0;i<3;i++){
        cells[count4][i] = cellfaces[i];
      }
       
      // printf("(%d, %d, %d, %d)\n", count4,cells[count4][0], cells[count4][1], cells[count4][2]);
      count2 += 1;
      count3 = 0;
      
    }
    
  }
  return(0);
}

int cmp(const void  *a, const void *b){
    int *x = *(int **)a;
    int *y = *(int **)b;
    // printf("%d %d %d\n",*x,*y, (*x > *y) - (*x < *y));
    if ((*x > *y) - (*x < *y) == 0){
        // printf("%d %d %d\n",*(x+1),*(y+1),(*(x+1) > *(y+1)) - (*(x+1) < *(y+1)));
        return (*(x+1) > *(y+1)) - (*(x+1) < *(y+1));
    }
    return (*x > *y) - (*x < *y);
}

int edge_comparison(int** edge_to_cell,int** edges,int** boundary_edges,int** cells,int** node_list,int npc,int* num_edges,int* num_boundary_edges,int skip,int edge_count){
  int i,j,count1,count2;
  count1 = 0;
  count2 = 0;
  for (i=0;i<edge_count;i++){       
    
    
    if (skip == 1){
      for (j=0;j<5;j++){
          node_list[0][j] = edge_to_cell[i][j];
          
      }
      
      
      skip = 0;
      continue;
    }

    for (j=0;j<5;j++){
        node_list[1][j] = edge_to_cell[i][j];
        
    }
    
    if (node_list[0][0]==node_list[1][0] && node_list[0][1]==node_list[1][1]){
    
        edges[count1][0] = node_list[0][0];
        edges[count1][1] = node_list[0][1];
        edges[count1][2] = node_list[0][2];
        edges[count1][3] = node_list[1][2];
        edges[count1][4] = 1;
        edges[count1][5] = 0;
        edges[count1][6] = node_list[0][3];
        edges[count1][7] = node_list[1][3];
        edges[count1][8] = 0;
        edges[count1][9] = node_list[0][4];
        edges[count1][10] = node_list[1][4];
        cells[node_list[0][2]][npc+node_list[0][3]] = count1;
        cells[node_list[1][2]][npc+node_list[1][3]] = count1;
        cells[node_list[0][2]][npc*2+node_list[0][3]] = node_list[1][2];
        cells[node_list[1][2]][npc*2+node_list[1][3]] = node_list[0][2];
        cells[node_list[0][2]][npc*3] = 0;
        cells[node_list[1][2]][npc*3] = 1;
        skip = 1;
        
        count1 += 1;
    }else{
          
        edges[count1][0] = node_list[0][0];
        edges[count1][1] = node_list[0][1];
        edges[count1][2] = node_list[0][2];
        edges[count1][3] = 0;
        edges[count1][4] = 5;
        edges[count1][5] = 1;
        edges[count1][6] = node_list[0][3];
        edges[count1][7] = 0;
        edges[count1][8] = count2;
        edges[count1][9] = node_list[0][4];
        edges[count1][10] = 0;
        for (j=0;j<11;j++){
            boundary_edges[count2][j] = edges[count1][j];
        }
        cells[node_list[0][2]][npc+node_list[0][3]] = count1;
        cells[node_list[0][2]][npc*2+node_list[0][3]] = 0;
        for (j=0;j<5;j++){
            node_list[0][j] = node_list[1][j];
        }
        
        skip = 0;
        count1 += 1;
        count2 += 1;
    }
  }
  if (skip == 0){

    edges[count1][0] = node_list[1][0];
    edges[count1][1] = node_list[1][1];
    edges[count1][2] = node_list[1][2];
    edges[count1][3] = 0;
    edges[count1][4] = 5;
    edges[count1][5] = 1;
    edges[count1][6] = node_list[1][3];
    edges[count1][7] = 0;
    edges[count1][8] = count2;
    edges[count1][9] = node_list[1][4];
    edges[count1][10] = 0;
    for (j=0;j<11;j++){
        boundary_edges[count2][j] = edges[count1][j];
    }
    cells[node_list[1][2]][npc+node_list[1][3]] = count1;
    cells[node_list[1][2]][npc*2+node_list[1][3]] = 0;
    count1 += 1;
    count2 += 1;
  }
  
  *num_edges = count1;
  *num_boundary_edges = count2;
  printf("Num edges: %d %d \n",*num_edges,*num_boundary_edges);
  return(0);
}

void freeMatrix_i2(int **matrix ,int row)
{
  for(int i=0;i<row;i++)
  {
    free(matrix[i]);
  }
  free(matrix);
}

void freeMatrix_d2(double **matrix ,int row)
{
  for(int i=0;i<row;i++)
  {
    free(matrix[i]);
  }
  free(matrix);
}

void freeMatrix_d3(double ***matrix ,int row,int col)
{
  for(int i=0;i<row;i++)
  {
    for(int j=0;j<col;j++){
      free(matrix[i][j]);
    }
    free(matrix[i]);
  }
  free(matrix);
}



double *** allocateMatrix_d3(int row,int col,int elems)
{
  double ***matrix = malloc(sizeof(double **)*row);
  for (int i = 0;i < row;i++){
    matrix[i] = malloc(sizeof(double*)*col);
    for (int j = 0;j < col;j++){
      matrix[i][j] = malloc(sizeof(double)*elems);
    }
  }
  return matrix;
}

int ** allocateMatrix_i2(int row,int elems)
{
  int **matrix = malloc(sizeof(int *)*row);
  for (int i = 0;i < row;i++){
    matrix[i] = malloc(sizeof(int)*elems);
  }
  return matrix;
}

double ** allocateMatrix_d2(int row,int elems)
{
  double **matrix = malloc(sizeof(double *)*row);
  for (int i = 0;i < row;i++){
    matrix[i] = malloc(sizeof(double)*elems);
  }
  return matrix;
}

double ** allocateMatrix_d2_zeros(int row,int elems)
{
  double **matrix = calloc(row,sizeof(double *));
  for (int i = 0;i < row;i++){
    matrix[i] = calloc(elems,sizeof(double));
  }
  return matrix;
}

void edge_quantities(int** edges,double** nodes,double** d_cell_face,double* area_face,double* area_bface,double** cell_center,double*** n̂_face,double*** face_mid,double** weight_face,double* weight_bface,double*** lf,double** lbf,double* δ_face,double* δ_bface,int num_edges,int vec_size)
{
  int ie,j,k;
  int node1,node2,cell1,cell2,face_cell1,face_cell2,bface_number;
  double scalar_result;
  double vector_result[3];
  for (ie=0;ie<num_edges;ie++){
    node1 = edges[ie][0];
    node2 = edges[ie][1];
    cell1 = edges[ie][2];
    cell2 = edges[ie][3];
    face_cell1 = edges[ie][6];
    face_cell2 = edges[ie][7];
    bface_number = edges[ie][8];
    vector_add_3(nodes[node2],1.0,nodes[node1],-1.0,vector_result);
    area_face[ie] = norm_3(vector_result);
    if (edges[ie][5] == 0){
      weight_face[ie][0] = d_cell_face[cell2][face_cell2]/(d_cell_face[cell1][face_cell1] + d_cell_face[cell2][face_cell2]);
      weight_face[ie][1] = 1.0 - weight_face[ie][0];

      vector_add_3(cell_center[cell2],1.0,cell_center[cell1],-1.0,lf[ie][0]);
        
      vector_mult_3(lf[ie][0], -1.0,lf[ie][1]);
        
      δ_face[ie] = fabs(dot_3(lf[ie][0],n̂_face[cell1][face_cell1]));     
    }
    else{
      area_bface[bface_number] = area_face[ie];
      weight_face[ie][0] = 1.0;
      weight_face[ie][1] = 0.0;
      weight_bface[bface_number] = 1.0;
      
      vector_add_3(face_mid[cell1][face_cell1],1.0,cell_center[cell1],-1.0,lf[ie][0]);
      vector_add_3(face_mid[cell1][face_cell1],1.0,cell_center[cell1],-1.0,lbf[bface_number]);
      δ_face[ie] = fabs(dot_3(lf[ie][0],n̂_face[cell1][face_cell1]));
      δ_bface[bface_number] = δ_face[ie];

    }
  

  }
  
}
void cell_quantities(int** circle_path,double*** vertices,int** cells,double** nodes,double*** t,double*** t̂,double*** face_mid,double** n_cell,double** n̂_cell,double*** n̂_face,double* volume_cell,double** d_cell_face,double** cell_center,int num_cells,int npc,int vec_size)
{
  int ic,j,k;
  double scalar_result;
  double vector_result[3];
  for (ic=0;ic<num_cells;ic++){
    for (j=0;j<npc;j++){
      for (k=0;k<vec_size;k++){
        vertices[ic][j][k] = nodes[cells[ic][j]][k];
        cell_center[ic][k] += vertices[ic][j][k];
      }
    }
    
    for (k=0;k<vec_size;k++){
      cell_center[ic][k] = cell_center[ic][k]/(double)npc;
      
    }
    
    for (j=0;j<npc;j++){
      for (k=0;k<vec_size;k++){
        t[ic][j][k] = vertices[ic][circle_path[j][1]][k] - vertices[ic][circle_path[j][0]][k];
        face_mid[ic][j][k] = (vertices[ic][circle_path[j][0]][k] + vertices[ic][circle_path[j][1]][k])/2.0;
      }
    }
    
    cross_product_3(t[ic][0],1.0,t[ic][1],-1.0,n_cell[ic]);
    scalar_result = norm_3(n_cell[ic]);
    volume_cell[ic] = scalar_result/2;
    
    for (k=0;k<vec_size;k++){
      n̂_cell[ic][k] = n_cell[ic][k]/scalar_result;
      
    }
    

      
    for (j=0;j<npc;j++){
      
      scalar_result = norm_3(t[ic][j]);
      for (k=0;k<vec_size;k++){
        t̂[ic][j][k] = t[ic][j][k]/scalar_result;
        
        
      }
      cross_product_3(t̂[ic][j],1.0,n̂_cell[ic],1.0,n̂_face[ic][j]);
      vector_add_3(cell_center[ic], 1.0, face_mid[ic][j], -1.0,vector_result);
      d_cell_face[ic][j] = norm_3(vector_result);
    }
    
    

  }
  

}

void node_weights(double** weight_node,double* weight_node_sum,double** nodes,int** cells,double** cell_center,int num_cells,int npc)
{
  double vector_result[3];
  for (int ic=0;ic<num_cells;ic++){
    for (int j=0;j<npc;j++){
      vector_add_3(nodes[cells[ic][j]],1.0,cell_center[ic],-1.0,vector_result);
      weight_node[ic][j] = 1.0/norm_3(vector_result);
      weight_node_sum[cells[ic][j]] += weight_node[ic][j];
    }
  }
  for (int ic=0;ic<num_cells;ic++){
    for (int j=0;j<npc;j++){
      weight_node[ic][j] = weight_node[ic][j]/weight_node_sum[cells[ic][j]];
    }
  }
}

void cross_product_3(double* a, double α_a, double* b, double α_b, double* c) {
  c[0] = α_a*a[1] * α_b*b[2] - α_a*a[2] * α_b*b[1];
  c[1] = α_a*a[2] * α_b*b[0] - α_a*a[0] * α_b*b[2];
  c[2] = α_a*a[0] * α_b*b[1] - α_a*a[1] * α_b*b[0];
}

void vector_add_3(double* a, double α_a, double* b, double α_b, double* c){
  c[0] = α_a*a[0] + α_b*b[0];
  c[1] = α_a*a[1] + α_b*b[1];
  c[2] = α_a*a[2] + α_b*b[2];
}

void vector_mult_3(double* a, double α_a, double* c){
  c[0] = α_a*a[0];
  c[1] = α_a*a[1];
  c[2] = α_a*a[2];
}

double norm_3(double* a){
    return sqrt(dot_3(a,a));

}

double dot_3(double* a,double* b){
  double c=0;
  for (int k=0;k<3;k++){
    c += a[k]*b[k];
  }
  return c;
}

void IC(int** cells,double* phi,double ICvalue,int num_cells){
    for (int ic=0;ic<num_cells;ic++){
      phi[ic] = ICvalue;
    }
}

void BC(int** boundary_edges,double* phi_bface,double* phi_node,struct BCarray* BCzone,int num_boundary_edges,int num_boundary_conditions){
  for (int ib=0;ib<num_boundary_edges;ib++){
    for (int bz=0;bz<num_boundary_conditions;bz++){
      if (boundary_edges[ib][4] == BCzone[bz].zone){
        if (BCzone[bz].type == 1){
          phi_bface[ib] = BCzone[bz].value;
        }
        else{
          if (BCzone[bz].type == 2){
            phi_bface[ib] = (phi_node[boundary_edges[ib][0]] + phi_node[boundary_edges[ib][1]])/2.0;
          }
        }
      }
    }
  }

}


void grid_properties(int** cells, int** edges,double* aₚ,double** aₙ,double* sourceᵤ,double* sourceₚ,double* sourceₛ,double** t̂_dot_lf,double*** t̂,double*** lf,double* area_face,double* δ_face,double Gamma,double Sourceᵤ,double Sourceₚ,double* volume_cell, struct BCarray* BCzone,int num_cells,int npc,int num_boundary_conditions)
{
  int ic,j,ibc;
  int edge;
  for (ic=0;ic<num_cells;ic++){
    aₚ[ic] = 0.0;
    sourceᵤ[ic] = Sourceᵤ * volume_cell[ic];
    sourceₚ[ic] = Sourceₚ * volume_cell[ic];
    for (j=0;j<npc;j++){
      edge = cells[ic][npc+j];
      if (edges[edge][5] == 0){
        aₙ[ic][j] = Gamma*area_face[edge]/δ_face[edge];
        aₚ[ic] += aₙ[ic][j];
        t̂_dot_lf[ic][j] = dot_3(t̂[ic][j],lf[edge][cells[ic][npc*3]]);
      }
      else{
        for (ibc=0;ibc<num_boundary_conditions;ibc++){
          aₙ[ic][j] = 0.0;
          t̂_dot_lf[ic][j] = dot_3(t̂[ic][j],lf[edge][0]);
          if (edges[edge][4] == BCzone[ibc].zone){
              if (BCzone[ibc].type == 2 ){
                  sourceᵤ[ic] += BCzone[ibc].value*area_face[edge];
              }else{
                  sourceᵤ[ic] += BCzone[ibc].value*Gamma*area_face[edge]/δ_face[edge];
                  sourceₚ[ic] += -Gamma*area_face[edge]/δ_face[edge];
              }
          }
        }
      }
    }
    aₚ[ic] -= sourceₚ[ic];
    
  }

  
  
}

void gauss_seidel(int** circle_path,int iterations,int** cells,double* phi_node,double* phi,double** weight_node,double* δ_face,double Gamma,int** boundary_edges,struct BCarray* BCzone,double* aₚ,double** aₙ,double* sourceᵤ,double* sourceₚ,double* sourceₛ,double** t̂_dot_lf,int npc,int num_cells,int num_boundary_edges,int num_boundary_conditions){
  int iter,ic,j,ib,ibc;
  for (iter = 0;iter<iterations;iter++){
    for (ic=0;ic<num_cells;ic++){
      for (j=0;j<npc;j++){
        phi_node[cells[ic][j]] += phi[ic]*weight_node[ic][j];
      }
    }
    for (ib=0;ib<num_boundary_edges;ib++){
      for (ibc=0;ibc<num_boundary_conditions;ibc++){
        if (boundary_edges[ib][4] == BCzone[ibc].zone){
          if (BCzone[ibc].type == 1 ){
            phi_node[boundary_edges[ib][0]] = BCzone[ibc].value;
            phi_node[boundary_edges[ib][1]] = BCzone[ibc].value;
          }
        } 
      }
    }
    for (ic=0;ic<num_cells;ic++){
      sourceₛ[ic] = 0.0;
      for (j=0;j<npc;j++){
          sourceₛ[ic] += t̂_dot_lf[ic][j]*Gamma/δ_face[cells[ic][j+npc]]*(phi_node[cells[ic][circle_path[j][1]]] - phi_node[cells[ic][circle_path[j][0]]]);
      }
    }
    for (ic=0;ic<num_cells;ic++){
      phi[ic] = (sourceᵤ[ic] + sourceₛ[ic])/aₚ[ic];
      for (j=0;j<npc;j++){
          if (cells[ic][j+npc*2] > 0){
              phi[ic] += aₙ[ic][j]*phi[cells[ic][j+npc*2]]/aₚ[ic];
              
          }
          
          phi_node[cells[ic][j]] = 0.0;
      }
    }
  }

}

void swap(int* arr, int i, int j) 
{ 
    int temp = arr[i]; 
    arr[i] = arr[j]; 
    arr[j] = temp; 
} 
  
// A function to implement bubble sort 
void bubbleSort(int arr[], int n) 
{ 
    int i, j; 
    for (i = 0; i < n - 1; i++) 
  
        // Last i elements are already 
        // in place 
        for (j = 0; j < n - i - 1; j++) 
            if (arr[j] > arr[j + 1]) 
                swap(arr, j, j + 1); 
} 