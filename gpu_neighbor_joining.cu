#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <climits>
#include <chrono>

using namespace std;
using namespace std::chrono;
#define TILE_WIDTH 16

const int MAX_TAXA = 100;

void checkCudaError(cudaError_t err)
{
    if (err != cudaSuccess) {
        printf("%s: %s\n", cudaGetErrorName(err), cudaGetErrorString(err));
        exit(1);
    }
}


struct node {
    int node_name;
    struct node* leftChild;
    struct node* rightChild;
    struct node* parent;
    double distance_left, distance_right;
};
typedef struct node Node;

Node* Node_new_all(int s, Node* lChild, Node* rChild, double lDistance, double rDistance){
    Node* this_node = (Node*)malloc(sizeof(Node));
    this_node->node_name     = s;
    this_node->leftChild     = lChild;
    this_node->rightChild    = rChild;
    this_node->parent        = nullptr;
    this_node->distance_left = lDistance;  
    this_node->distance_right = rDistance;  
    lChild->parent      = this_node;
    rChild->parent      = this_node;
    return this_node;
}

Node* Node_new(int s) {
    Node* this_node = (Node*)malloc(sizeof(Node));
    this_node->node_name     = s;
    this_node->leftChild     = nullptr;
    this_node->rightChild    = nullptr;
    this_node->parent        = nullptr;
    this_node->distance_left = -1;  
    this_node->distance_right = -1;  
    return(this_node);
}

int readFromFile(double* dist_mat, char seq[MAX_TAXA], string filename, Node* nodes[MAX_TAXA]) {
    cout<<filename<<endl;
    ifstream infile(filename);

    if (!infile) {
        cerr << "Error opening file" << endl;
        exit(1);
    }
    int num_taxa;
    infile >> num_taxa;
    infile.peek();
    int numRows = 0, numCols = 0;
    for (int i = 0; i < num_taxa; i++) {
        for (int j = 0; j < num_taxa; j++) {
            dist_mat[i*num_taxa + j] = 0;
        }
    }

    string s;
    int l_node = 0;
    while (!infile.eof() && numRows < num_taxa) {
        numCols = 0;
	    infile >> s;
        nodes[numRows] = Node_new(l_node);
	    infile.peek();
        while (infile.peek() != '\n' && numCols < numRows) {
            infile >> dist_mat[numRows*num_taxa + numCols];
            dist_mat[numCols*num_taxa + numRows] = dist_mat[numRows*num_taxa + numCols];
            numCols++;
        }
        infile.ignore(); // ignore newline character
        numRows++;
        l_node++;
    }

    infile.close();
    return num_taxa;
}

void printDistanceMatrix(double* dist_mat, int num_taxa, Node* nodes[MAX_TAXA]){
	cout<< "Num_taxa = " << num_taxa <<endl;
    for (int i = 0; i < num_taxa; i++) {
        if(nodes[i]==nullptr)        {
		    cout<<"Seq "<<i<<" = "<<"NULL" << " : ";
        }
        else {
		    cout<<"Seq "<<i<<" = "<<to_string(nodes[i]->node_name) << " : ";
        }
        for (int j = 0; j < num_taxa; j++) {
            cout << dist_mat[i*num_taxa + j] << " ";
        }
        cout << endl;
    }
}

void printTDMatrix(double TD_arr[MAX_TAXA], int num_taxa){
    for(int i=0 ; i<num_taxa; i++){
        cout<<TD_arr[i] << " " ;
    }
    cout<<endl;
}

///brief To dump values from Tree into a text file for graph generation
void traverseAndWrite(Node* node, ofstream& outfile) {
    if (node != NULL) {
        // Process the current node
        if(node->leftChild!=nullptr) {
            outfile << "\"" << to_string(node->node_name) << "\" ";
            outfile << "->" << "\"" << to_string(node->leftChild->node_name) << "\" ";
            outfile <<"[taillabel = " <<fixed << setprecision(2) << node->distance_left <<"]"<<endl;

        }
        if(node->rightChild!=nullptr) {
            outfile << "\"" << to_string(node->node_name) << "\" ";
            outfile << "->" << "\"" << to_string(node->rightChild->node_name) << "\" ";
            outfile <<"[taillabel = " <<fixed << setprecision(2) << node->distance_right <<"]"<<endl;

        }
        // Traverse the left subtree
        traverseAndWrite(node->leftChild, outfile);

        // Traverse the right subtree
        traverseAndWrite(node->rightChild, outfile);
    }
}



__global__ void gpu_nj(int num_taxa, double* d_dist_mat, double* d_TD_arr, double* d_TB_min, Node** d_nodes, Node* d_temp_node){

    /*
        Device variables needed - dist_mat, TD_arr
    */ 

    int tid = blockIdx.x*blockDim.x + threadIdx.x;
    int t_row = blockIdx.x*blockDim.x + threadIdx.y;
    int t_row_z = threadIdx.y;
    int t_col_og = threadIdx.x;
    int t_col = t_col_og;
    //int t_row = tid / num_taxa;
    //int t_col = tid % num_taxa;
    //__shared__ int index1, index2;
    //__shared__ int min_index, max_index;
    //__shared__ double delta_ij, limb_length_i, limb_length_j;
    __shared__ double D_star_mat[32][2]; 
    //__shared__ double s_td_arr[32];

    int index1, index2;
    __shared__ int min_index;
    __shared__ int max_index;
    __shared__ double delta_ij;
    __shared__ double limb_length_i;
    __shared__ double limb_length_j;
    //__shared__ double min_D_star_mat[TILE_WIDTH][2];
    __shared__ int col_active[TILE_WIDTH];
    __shared__ int row_active[TILE_WIDTH];
    __shared__ double sum_tile[TILE_WIDTH][TILE_WIDTH];
    //double d_TD_arr[32];
    __shared__ double min_row_mat[2][TILE_WIDTH][TILE_WIDTH];
    int n;
    int i, j;
    double min_d_star_row;



            
    int col_tile_iter;
    int par_sum_iter = TILE_WIDTH / 2;
   
    // FIXME: OPT - can go down the column per thread
    // parallel sum possible 
    for(i=0 ; i<num_taxa-2; i++) {
        n = num_taxa - i;

        // initialize sum to 0
        
        sum_tile[t_row_z][t_col_og] = 0;

        if(t_col_og == 0){
            row_active[t_row_z] = -1;
            if(d_dist_mat[t_row] != -1 && t_row < num_taxa){
                row_active[t_row_z] = 1;
            }
        }
        __syncthreads();
        
        // loop needed to iterate over the columns of dist_mat
        for(col_tile_iter = 0; col_tile_iter < gridDim.x; col_tile_iter++) {
            t_col = t_col_og + col_tile_iter*TILE_WIDTH;
            // load which col is active amongst the current TILE_WIDTH
            if(t_row_z == 0) {
                col_active[t_col_og] = -1;
                if(d_dist_mat[t_col] != -1 && t_col < num_taxa){
                    col_active[t_col_og] = 1;
                }
            }
            __syncthreads();
            // find the sum corresponding to current TILE_WIDTH
            if(col_active[t_col_og] == 1  && row_active[t_row_z] == 1){
                sum_tile[t_row_z][t_col_og] += d_dist_mat[t_row*num_taxa + t_col];
                //printf("row %d of sum_tile added\n", t_row_z);
            }
        }
        __syncthreads();

//if(t_row == 0 && t_col_og == 0) {
//    printf("Printing sum_tile_Arr\n");
//    for(int dbg = 0; dbg < TILE_WIDTH; dbg++){
//        for(int dbg2 = 0; dbg2 < TILE_WIDTH; dbg2++){
//            printf("%lf ", sum_tile[dbg][dbg2]);
//        }
//        printf("\n");
//    }
//    printf("\n");
//}

//if(t_row == 0 && t_col == 0) {
//    printf("Printing col_active\n");
//    for(int dbg = 0; dbg < TILE_WIDTH; dbg++){
//        printf("%d ", col_active[dbg]);
//    }
//    printf("\n");
//}
//if(t_row == 0 && t_col == 0) {
//    printf("Printing row_active\n");
//    for(int dbg = 0; dbg < TILE_WIDTH; dbg++){
//        printf("%d ", row_active[dbg]);
//    }
//    printf("\n");
//}
//__syncthreads();

        
        // now find the sum of that row of the memory with parallel sum
        /*
        while (par_sum_iter > 0) {
            if(t_col_og < par_sum_iter){
                sum_tile[t_row_z][t_col_og] += sum_tile[t_row_z][t_col_og + par_sum_iter];
                //printf("par_sum_iter = %d\n", par_sum_iter);
            }
            __syncthreads();
            par_sum_iter = par_sum_iter / 2;
        }
        
        if(t_col_og == 0){
            if(col_active[t_row_z] == 1){
                d_TD_arr[t_row] = sum_tile[t_row_z][0];
            }
            else {
                d_TD_arr[t_row] = -1;
            }
        }
        */
        
    
        // store the calculated sum to d_TD_arr
        
        if(t_col_og == 0){
            for(par_sum_iter = 1; par_sum_iter < TILE_WIDTH; par_sum_iter++){
                sum_tile[t_row_z][0] += sum_tile[t_row_z][t_col_og + par_sum_iter];
            }
            if(row_active[t_row_z] == 1){
                d_TD_arr[t_row] = sum_tile[t_row_z][0];
            }
            else {
                d_TD_arr[t_row] = -1;
            }
        }
        
//        }
        __syncthreads();

        // DEBUG print the TD_arr

//if(t_row == 0 && t_col == 0) {
//    printf("Printing TD_Arr\n");
//    for(int dbg = 0; dbg < num_taxa; dbg++){
//        printf("%lf ", d_TD_arr[dbg]);
//    }
//    printf("\n");
//}
//__syncthreads();        




        
        // initialize min_row_mat with INT_MAX
        min_row_mat[0][t_row_z][t_col_og] = INT_MAX;
        min_row_mat[1][t_row_z][t_col_og] = t_col_og;

        double curr_value, calc_value, curr_index;
        curr_value = INT_MAX;
        curr_index = t_col_og;
        __syncthreads();

        // for each tile iterate and find the min of the D_star similar to sum
        
        for(col_tile_iter = 0; col_tile_iter < gridDim.x; col_tile_iter++) {
            t_col = t_col_og + col_tile_iter*TILE_WIDTH;
            // load which col is active amongst the current TILE_WIDTH
            if(t_row_z == 0) {
                col_active[t_col_og] = -1;
                if(d_dist_mat[t_col] != -1 && t_col < num_taxa){
                    col_active[t_col_og] = 1;
                }
            }
            __syncthreads();
            //  find the min of current and iter 
            if(col_active[t_col_og] == 1 && row_active[t_row_z] == 1){
                // calculate the D_star value for the elem
                // compare with the existing min
                // replace if less
                calc_value = (n - 2) * d_dist_mat[t_row*num_taxa + t_col] - d_TD_arr[t_row] - d_TD_arr[t_col];
                if(calc_value < curr_value){
                    curr_value = calc_value;
                    curr_index = t_col;
                }
            }
        }
        // store the curr_index and curr_value
        min_row_mat[0][t_row_z][t_col_og] = curr_value;
        min_row_mat[1][t_row_z][t_col_og] = curr_index;
        __syncthreads();
        
        int min_iter;
        // find the min of the row and its index
        if(t_col_og == 0 && row_active[t_row_z] == 1){
            for(min_iter = 1; min_iter < TILE_WIDTH; min_iter++){
                if(min_row_mat[0][t_row_z][min_iter] < curr_value) {
                   curr_value = min_row_mat[0][t_row_z][min_iter];
                   curr_index = min_row_mat[1][t_row_z][min_iter];
                }
            }
            min_row_mat[0][t_row_z][0] = curr_value;
            min_row_mat[1][t_row_z][0] = curr_index;
        }
        __syncthreads();
        // PRINT THE MIN_ROW_MAT[0] values                

        // find the min of all rows and the index pair
        if(t_row_z == 0 && t_col_og == 0){
            for(min_iter = 0; min_iter < TILE_WIDTH; min_iter++){
                if(min_row_mat[0][min_iter][0] < curr_value){
                    curr_value = min_row_mat[0][min_iter][0];
                    curr_index = min_row_mat[1][min_iter][1];
                }
            }
            d_TB_min[blockIdx.x] = curr_value;
            d_TB_min[gridDim.x + blockIdx.x] = t_row; // row
            d_TB_min[2*gridDim.x + blockIdx.x] = curr_index; // col
        }
        __syncthreads();
                        
        // reduce the d_TB_min
        if(t_row == 0 && t_col_og == 0){
            curr_value = d_TB_min[0];
            index1 = d_TB_min[gridDim.x];
            index2 = d_TB_min[2*gridDim.x];
            for(min_iter = 1; min_iter < gridDim.x; min_iter++){
                if(d_TB_min[min_iter] < curr_value){
                    curr_value = d_TB_min[min_iter];
                    index1 = d_TB_min[gridDim.x + min_iter];
                    index2 = d_TB_min[2*gridDim.x + min_iter];
                }
            }
        }
                                     
                                     
                                     
                                     
                                     
                                     
        /*                           
        //find_closest_pair(d        _dist_mat,num_taxa, TD_arr, index1, index2);
        // GPU code for find_        closest_pair
        min_d_star_row = INT_        MAX;
                                     
        if(tid < num_taxa) {         
            D_star_mat[tid][0        ] = INT_MAX;
            if(d_dist_mat[tid        ] != -1) {
                for (j = tid         + 1; j < num_taxa; j++) {
                    if(d_dist        _mat[j] != -1){
                        min_d        _star_row = (n - 2) * d_dist_mat[tid*num_taxa + j] - d_TD_arr[tid] - d_TD_arr[j];
                        if (min_d_star_row < D_star_mat[tid][0]) {
                            D_star_mat[tid][0] = min_d_star_row;
                            D_star_mat[tid][1] = j;
                            
                            }
                        }
                    }
                } else {
                    D_star_mat[tid][0] = INT_MAX;
                }
            }
        __syncthreads();
        

        // find the index pair which has absolute min among the d_star
        if(tid == 0) {
            min_d_star_row = INT_MAX;
            for (j = 0; j < num_taxa; j++) {
                if(D_star_mat[j][0] < min_d_star_row){
                    min_d_star_row = D_star_mat[j][0];
                    index1 = j;
                    index2 = D_star_mat[j][1];
                }
            }
        }
        */




        if(tid == 0) {
            min_index = (index1 < index2) ? index1 : index2;
            max_index = (index1 < index2) ? index2 : index1;
            delta_ij = (d_TD_arr[min_index] - d_TD_arr[max_index]) / (n-2);
            limb_length_i = (d_dist_mat[min_index*num_taxa + max_index] + delta_ij)/2.0;
            limb_length_j = (d_dist_mat[min_index*num_taxa + max_index] - delta_ij)/2.0;
        }

        __syncthreads();


        //updateDistanceMatrix(d_dist_mat,num_taxa, min_index, max_index);

        // update the distance matrix parallely with new values at max index

        if((tid < num_taxa) && (tid != min_index) && (tid != max_index)) {
            d_dist_mat[max_index*num_taxa + tid] = (d_dist_mat[min_index*num_taxa + tid] + d_dist_mat[max_index*num_taxa + tid] - d_dist_mat[min_index*num_taxa + max_index]) / 2;
            //printf("Modified value of d_dist_mat[%d][%d] to %lf \n", max_index, tid, d_dist_mat[max_index*num_taxa + tid]);
            d_dist_mat[tid*num_taxa + max_index] = d_dist_mat[max_index*num_taxa + tid];
        }


        // turn min_index to -1,
        if(tid == 0){
            d_dist_mat[min_index*num_taxa] = -1;
            d_dist_mat[min_index] = -1;
         
            d_temp_node->node_name      = i;
            d_temp_node->leftChild      = d_nodes[min_index];
            d_temp_node->rightChild     = d_nodes[max_index];
            d_temp_node->parent         = nullptr;
            d_temp_node->distance_left  = limb_length_i;
            d_temp_node->distance_right = limb_length_j;

            // should not create a new node in GPU, rather just change the values of existing array
            d_nodes[max_index] = d_temp_node;
            d_nodes[min_index] = nullptr;
        }
        __syncthreads();
        
        //print the matrix
        /*
        if(tid == 0){
            for(int p=0; p<num_taxa; p++) {
                printf("%lf ",d_TD_arr[p]);
            }
            
            printf(" \n");
            printf("min_d_star is %lf for i=%d, j=%d\n", min_d_star_row, min_index, max_index);
            printf("delta_ij = %lf, Li = %lf, Lj = %lf\n", delta_ij, limb_length_i, limb_length_j);

            printf("-- Printing d_dist_mat for iteration %d -- \n", i);
            for(int p=0; p<num_taxa; p++){
                printf(" row %d ", p);
                for(int q=0; q<num_taxa; q++) {
                    printf("%lf ", d_dist_mat[p*num_taxa + q]);
                }
                printf(" \n");
            }

        }
        */
    }        
             
    // Copyi ng the TD_arr back to GPU global memory

    if(tid < num_taxa);
        d_TD_arr[tid] = d_TD_arr[tid];

    //__syncthreads();

};



int main() {
    
    //string filename = "./examples/evolution.in";
    string filename = "./examples/INGI2368.in";
    ifstream infile(filename);
    if (!infile) {
        cerr << "Error opening file" << endl;
        exit(1);
    }
    int num_taxa;
    infile >> num_taxa;
    //double dist_mat[num_taxa*num_taxa];
    double* dist_mat;
    dist_mat = (double *)malloc(num_taxa*num_taxa*sizeof(double));
    double* d_dist_mat;
    char seq[num_taxa];
    Node* nodes[num_taxa];
    readFromFile(dist_mat, seq, filename, nodes);
    Node** d_nodes;
    Node* d_temp_node;
    printDistanceMatrix(dist_mat, num_taxa, nodes);
    //int index1, index2;
    //int min_index, max_index;
    //double delta_ij, limb_length_i, limb_length_j;
    //int n;
    double TD_arr[num_taxa];
    double* d_TD_arr;
    double* d_TB_min;
    int num_TB = (num_taxa + TILE_WIDTH - 1) / TILE_WIDTH;

    // allocate memory and copy the variables to GPU, 
    // launch kernel
    // copy the variables to CPU
    // free GPU memory

    printf("*** Allocating GPU memory ***\n");
    cudaMalloc((void**)(&d_dist_mat), num_taxa*num_taxa*(sizeof(double)));
    cudaMalloc((void**)(&d_TD_arr), num_taxa*(sizeof(double)));
    cudaMalloc((void**)(&d_TB_min), 3*num_TB*(sizeof(double)));
    cudaMalloc((void**)(&d_nodes), num_taxa*(sizeof(Node)));
    cudaMalloc((void**)(&d_temp_node), sizeof(Node));
    printf("*** Allocating GPU memory complete ***\n\n");

    printf("*** Copying to GPU memory ***\n");
    checkCudaError(cudaMemcpy(d_dist_mat, dist_mat, num_taxa*num_taxa*(sizeof(double)), cudaMemcpyHostToDevice));    
    //cudaMemcpy(&d_TD_arr, &TD_arr, num_taxa*(sizeof(double)), cudaMemcpyHostToDevice);
    //cudaMemcpy(&d_TB_min, &d_TB_min, 3*num_TB*(sizeof(double)), cudaMemcpyHostToDevice);
    //cudaMemcpy(&d_nodes, &nodes, num_taxa*(sizeof(Node)), cudaMemcpyHostToDevice);
    printf("*** Copying to GPU memory complete ***\n\n");

    // Parallelize GPU set grid, block and call kernel
    dim3 blocksize(TILE_WIDTH,TILE_WIDTH);
    dim3 gridsize((num_taxa + TILE_WIDTH - 1) / TILE_WIDTH);
   
    //printf("Launching kernel with griddim: %d, %d, %d\n", gridDim.x, gridDim.y, gridDim.z);
    //printf("Launching kernel with blockDim: %d, %d, %d\n", blockDim.x, blockDim.y, blockDim.z);
    
    auto start_time = high_resolution_clock::now(); 
    gpu_nj<<<gridsize, blocksize>>>(num_taxa, d_dist_mat, d_TD_arr, d_TB_min,  d_nodes, d_temp_node);
    cudaDeviceSynchronize();    
    auto end_time = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end_time - start_time);
    printf("### \n Elapsed Time %" PRId64 "\n###\n", duration.count());
   
 
    //checkCudaError(cudaGetLastError());
    printf("***  GPU computation complete ***\n");

    checkCudaError(cudaMemcpy(dist_mat, d_dist_mat, num_taxa*num_taxa*sizeof(double), cudaMemcpyDeviceToHost));
    printf("*** Transferring data from Device to Host complete ***\n");
    checkCudaError(cudaMemcpy(TD_arr, d_TD_arr, num_taxa*sizeof(double), cudaMemcpyDeviceToHost));
    printf("*** Transferring data from Device to Host complete ***\n");
    //checkCudaError(cudaMemcpy(nodes, *d_nodes, num_taxa*sizeof(Node), cudaMemcpyDeviceToHost));
    printf("*** Transferring data from Device to Host complete ***\n");

    int final_index1 = -1;
    int final_index2 = -1;

    int i;
    for(i=0 ; i<num_taxa ; i++) {
        if(dist_mat[i*num_taxa + 0]!=-1)
        {
            if(final_index1==-1)
                final_index1 = i;
            else
                final_index2 = i;
        }
    } 

    int root_node_name = i;
    cout<<to_string(root_node_name)<<endl;
    Node* root = Node_new_all(root_node_name, nodes[final_index1], nodes[final_index2], dist_mat[final_index1*num_taxa + final_index2]/2.0, dist_mat[final_index1*num_taxa + final_index2]/2.0 );
    
    printf("*** Final node computed ***\n");
    printDistanceMatrix(dist_mat, num_taxa, nodes);

    ofstream outfile("g.gv"); // open the output file
    if (!outfile) {
        cerr << "Error opening file" << endl;
        exit(1);
    }
    outfile << "digraph {" << endl;
    traverseAndWrite(root, outfile);
    outfile << "}" << endl;
    outfile.close();

    //FIXME:free the gpu memory and all variables

    cudaFree(d_dist_mat);
    cudaFree(d_TD_arr);
    cudaFree(d_TB_min);
    cudaFree(d_nodes);
    cudaFree(d_temp_node);

    return 0;
}
