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

using namespace std;

const int MAX_TAXA = 100;

struct node {
    public:
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

int readFromFile(double dist_mat[MAX_TAXA][MAX_TAXA], char seq[MAX_TAXA], string filename, Node* nodes[MAX_TAXA]) {
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
    //Initialize Distance Matrix to 0
    for (int i = 0; i < num_taxa; i++) {
        for (int j = 0; j < num_taxa; j++) {
            dist_mat[i][j] = 0;
        }
    }

    while (!infile.eof() && numRows < num_taxa) {
        numCols = 0;
	    infile >> seq[numRows];
        nodes[numRows] = Node_new({seq[numRows]});
	    infile.peek();
        while (infile.peek() != '\n' && numCols < numRows) {
            infile >> dist_mat[numRows][numCols];
            dist_mat[numCols][numRows] = dist_mat[numRows][numCols];
            numCols++;
        }
        infile.ignore(); // ignore newline character
        numRows++;
    }

    infile.close();
    return num_taxa;
}

void printDistanceMatrix(double dist_mat[MAX_TAXA][MAX_TAXA], int num_taxa, Node* nodes[MAX_TAXA]){
	cout<< "Num_taxa = " << num_taxa <<endl;
    for (int i = 0; i < num_taxa; i++) {
        if(nodes[i]==nullptr)        {
		    cout<<"Seq "<<i<<" = "<<"NULL" << " : ";
        }
        else {
		    cout<<"Seq "<<i<<" = "<<to_string(nodes[i]->node_name) << " : ";
        }
        for (int j = 0; j < num_taxa; j++) {
            cout << dist_mat[i][j] << " ";
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


void totalDistance(double dist_mat[MAX_TAXA][MAX_TAXA], int num_taxa, double TD_arr[MAX_TAXA]){

    // total NUM_TAXA threads,
    // each thread will be calculating the sum for TD_arr[i]
    for(int i=0; i<num_taxa; i++){
        double sum=0;
        if(dist_mat[i][0]!=-1) {
            for (int k = 0; k < num_taxa; k++) {
                if(dist_mat[k][0]!=-1){
                    sum += dist_mat[i][k];
                }
            }
            TD_arr[i] = sum;
        } else{
            TD_arr[i] = -1;
        }
    }
}

/// @brief Calculate indexes with minimum D_star and store in array variable pair
// Check if -1
void find_closest_pair(double dist_mat[MAX_TAXA][MAX_TAXA], int num_taxa, double TD_arr[MAX_TAXA], int& index1, int& index2) {
    // less than num_taxa*num_taxa/2 threads used
    // change this code:
    // first create a D_star matrix in shared memory - each thread independently
    // Then find the minimum of the D_star again using parallelism using parallel_red
    
    double min_distance = INT_MAX;
    for (int i = 0; i < num_taxa; i++) {
        if(dist_mat[i][0]!=-1) {
            for (int j = i + 1; j < num_taxa; j++) {
                if(dist_mat[j][0]!=-1){
                    double D_star = (num_taxa - 2) * dist_mat[i][j] - TD_arr[i] - TD_arr[j];
                    if (D_star < min_distance) {
                        min_distance = D_star;
                        index1 = i;
                        index2 = j;
                    }
                }
            }
        }
    }
    //return min_distance / (num_taxa - 2);
}

void updateDistanceMatrix(double dist_mat[MAX_TAXA][MAX_TAXA], int num_taxa, int min_index, int max_index) {
    // update the distance matrix parallely with new values
    for (int k = 0; k < num_taxa; k++) {
        if (k != min_index && k != max_index) {
            dist_mat[max_index][k] = ( dist_mat[min_index][k] + dist_mat[max_index][k] - dist_mat[min_index][max_index]) / 2;
            dist_mat[k][max_index] = dist_mat[max_index][k];
        }
    }
    dist_mat[min_index][0] = dist_mat[0][min_index] = -1;
}

__global__ void gpu_nj(num_taxa, dist_mat, TD_arr, temp_dist_mat){

    /*
        Device variables needed - dist_mat, TD_arr
    */ 

    int tid = blockIdx.x*blockDim.x + threadIdx.x;
    int t_row = tid / num_taxa;
    int t_col = tid % num_taxa;
    __shared__ int index1, index2;
    __shared__ int min_index, max_index;
    __shared__ double delta_ij, limb_length_i, limb_length_j;
    int n;
    int i;
    double sum, min_d_star_row;
    int new_node_name;
    Node* temp_node;

    __shared__ double D_star_mat[100][2]; // declare in shared memory later

    // load dist_mat in shared memory
    // load TD_arr in shared memory
    
   
    // OPT - can go down the column per thread
    // parallel sum possible 
    for(int i=0 ; i<num_taxa-2; i++) {
        n = num_taxa - i;
        //totalDistance(dist_mat, num_taxa, TD_arr);
        // GPU implementation of totalDistance
        if(tid < num_taxa) {
            if(dist_mat[tid][0] != -1) {
                sum=0;
                for (int k = 0; k < num_taxa; k++) {
                    if(dist_mat[k][0] != -1){
                        sum += dist_mat[tid][k];
                    }
                }
                TD_arr[tid] = sum;
            } else{
                TD_arr[tid] = -1;
            }
        }


        //find_closest_pair(dist_mat,num_taxa, TD_arr, index1, index2);
        // GPU code for find_closest_pair
        min_d_star_row = INT_MAX;
        if(tid < num_taxa) {
            if(dist_mat[tid][0] != -1) {
                for (int j = tid + 1; j < num_taxa; j++) {
                    if(dist_mat[j][0] != -1){
                        min_d_star_row = (num_taxa - 2) * dist_mat[tid][j] - TD_arr[tid] - TD_arr[j];
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
        }

        // find the index pair which has absolute min among the d_star
        if(tid == 0) {
            min_d_star_row = INT_MAX;
            for (i = 0; i < num_taxa; i++) {
                if(D_star_mat[tid][0] < min_d_star_row){
                    min_d_star_row = D_star_mat[tid][0];
                    index1 = i;
                    index2 = D_star_mat[tid][1];
                }
            }
        }


        if(tid == 0) {
            min_index = index1 < index2 ? index1 : index2;
            max_index = index1 < index2 ? index2 : index1;
            delta_ij = (TD_arr[min_index] - TD_arr[max_index]) / (n-2);
            limb_length_i = (dist_mat[min_index][max_index] + delta_ij)/2.0;
            limb_length_j = (dist_mat[min_index][max_index] - delta_ij)/2.0;
        }


        //updateDistanceMatrix(dist_mat,num_taxa, min_index, max_index);

        // update the distance matrix parallely with new values
       
        // create a new distance matrix and swap the pointers 
        if((t_row != min_index) && (t_col !=  min_index) && (t_row != max_index) && (t_col !=  max_index)) {
                temp_dist_mat[t_row][t_col] = (dist_mat[min_index][t_col] + dist_mat[max_index][t_col] - dist_mat[min_index][max_index]) / 2;
        temp_dist_mat[min_index][0] = dist_mat[0][min_index] = -1;
        }
        if(tid == 0){
            temp_dist_mat[min_index][0] = -1;
            temp_dist_mat[0][min_index] = -1;
            dist_mat = temp_dist_mat;
   
         
            new_node_name = i;
            temp_node = Node_new_all(new_node_name, nodes[min_index], nodes[max_index], limb_length_i, limb_length_j );
            nodes[max_index] = temp_node;
            nodes[min_index] = nullptr;
        }
        __syncthreads();
    }

};



int main() {
    
    string file_name = "./examples/INGI2368.in";
    double dist_mat[MAX_TAXA][MAX_TAXA];
    char seq[MAX_TAXA];
    Node* nodes[MAX_TAXA];
    //int num_taxa = read_DM_file(dist_mat, seq, file_name, nodes);
    int num_taxa = readFromFile(dist_mat, seq, file_name, nodes);
    printDistanceMatrix(dist_mat, num_taxa, nodes);
    int index1, index2;
    int min_index, max_index;
    double delta_ij, limb_length_i, limb_length_j;
    int n;
    double TD_arr[MAX_TAXA];

    // Parallelize GPU
    for(int i=0 ; i<num_taxa-2; i++) {
        n = num_taxa - i;
        totalDistance(dist_mat, num_taxa, TD_arr);
        printTDMatrix(TD_arr, num_taxa);
        find_closest_pair(dist_mat,num_taxa, TD_arr, index1, index2);
        
        min_index = min(index1, index2);
        max_index = max(index1, index2);
        delta_ij = (TD_arr[min_index] - TD_arr[max_index]) / (n-2);
        limb_length_i = (dist_mat[min_index][max_index] + delta_ij)/2.0;
        limb_length_j = (dist_mat[min_index][max_index] - delta_ij)/2.0;
        updateDistanceMatrix(dist_mat,num_taxa, min_index, max_index);
        int new_node_name = i;
        cout<<to_string(new_node_name)<<endl;
        Node* temp = Node_new_all(new_node_name, nodes[min_index], nodes[max_index], limb_length_i, limb_length_j );
        nodes[max_index] = temp;
        nodes[min_index] = nullptr;
        printDistanceMatrix(dist_mat, num_taxa, nodes);
    }



    int final_index1 = -1;
    int final_index2 = -1;

    int i;
    for(i=0 ; i<num_taxa ; i++) {
        if(dist_mat[i][0]!=-1)
        {
            if(final_index1==-1)
                final_index1 = i;
            else
                final_index2 = i;
        }
    } 

    int root_node_name = i;
    cout<<to_string(root_node_name)<<endl;
    Node* root = Node_new_all(root_node_name, nodes[final_index1], nodes[final_index2], dist_mat[final_index1][final_index2]/2.0, dist_mat[final_index1][final_index2]/2.0 );

    // cout<<nodes[0]->node_name<<" "<<nodes[1]->node_name;
    
    ofstream outfile("g.gv"); // open the output file
    if (!outfile) {
        cerr << "Error opening file" << endl;
        exit(1);
    }
    outfile << "digraph {" << endl;
    traverseAndWrite(root, outfile);
    outfile << "}" << endl;
    outfile.close();


    return 0;
}
