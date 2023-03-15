#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <climits>
#include <chrono>


using namespace std;
using namespace std::chrono;

const int MAX_TAXA = 100;

class Node {
    public:
    string node_name;
    Node* leftChild;
    Node* rightChild;
    Node* parent;
    double distance_left, distance_right;
    /**
     * Node constructor, which initializes everything
     */
    Node(string s) {
        this->node_name     = s;
        this->leftChild     = nullptr;
        this->rightChild    = nullptr;
        this->parent        = nullptr;
        this->distance_left = -1;  
        this->distance_right = -1;  
    }

    Node(string s, Node* lChild, Node* rChild, double lDistance, double rDistance) {
        this->node_name     = s;
        this->leftChild     = lChild;
        this->rightChild    = rChild;
        this->parent        = nullptr;
        this->distance_left = lDistance;  
        this->distance_right = rDistance;  
        lChild->parent      = this;
        rChild->parent      = this;
    }
};

int readFromFile(double arr[MAX_TAXA][MAX_TAXA], char seq[MAX_TAXA], string filename, Node* nodes[MAX_TAXA]) {
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
            arr[i][j] = 0;
        }
    }

    string s;
    while (!infile.eof() && numRows < num_taxa) {
        numCols = 0;
	    infile >> s;
        nodes[numRows] = new Node(s);
	    infile.peek();
        while (infile.peek() != '\n' && numCols < numRows) {
            infile >> arr[numRows][numCols];
            arr[numCols][numRows] = arr[numRows][numCols];
            numCols++;
        }
        infile.ignore(); // ignore newline character
        numRows++;
    }

    infile.close();
    return num_taxa;
}

void printDistanceMatrix(double arr[MAX_TAXA][MAX_TAXA], int num_taxa, Node* nodes[MAX_TAXA]){
	cout<< "Num_taxa = " << num_taxa <<endl;
    for (int i = 0; i < num_taxa; i++) {
        if(nodes[i]==nullptr)        {
		    cout<<"Seq "<<i<<" = "<<"NULL" << " : ";
        }
        else {
		    cout<<"Seq "<<i<<" = "<<nodes[i]->node_name << " : ";
        }
        for (int j = 0; j < num_taxa; j++) {
            cout << arr[i][j] << " ";
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
// void getRandomIndexes(int& random_number1, int& random_number2, int range ){
//     std::srand(std::time(nullptr));
//     random_number1 = std::rand() % range;
//     random_number2 = std::rand() % range;
//     if(random_number2==random_number1)
//         random_number2 = (random_number1 + 1) % range;
// }

///brief To dump values from Tree into a text file for graph generation
void traverseAndWrite(Node* node, ofstream& outfile) {
    if (node != NULL) {
        // Process the current node
        if(node->leftChild!=nullptr) {
            outfile << "\"" << node->node_name << "\" ";
            outfile << "->" << "\"" << node->leftChild->node_name << "\" ";
            outfile <<"[taillabel = " <<fixed << setprecision(2) << node->distance_left <<"]"<<endl;

        }
        if(node->rightChild!=nullptr) {
            outfile << "\"" << node->node_name << "\" ";
            outfile << "->" << "\"" << node->rightChild->node_name << "\" ";
            outfile <<"[taillabel = " <<fixed << setprecision(2) << node->distance_right <<"]"<<endl;

        }
        // Traverse the left subtree
        traverseAndWrite(node->leftChild, outfile);

        // Traverse the right subtree
        traverseAndWrite(node->rightChild, outfile);
    }
}


void totalDistance(double arr[MAX_TAXA][MAX_TAXA], int num_taxa, double TD_arr[MAX_TAXA]){
    for(int i=0; i<num_taxa; i++){
        double sum=0;
        TD_arr[i] = -1;
        if(arr[i][0]!=-1) {
            for (int k = 0; k < num_taxa; k++) {
                if(arr[k][0]!=-1){
                    sum += arr[i][k];
                }
            }
            TD_arr[i] = sum;
        }
    }
}

/// @brief Calculate indexes with minimum D_star and store in array variable pair
// Check if -1
void find_closest_pair(double arr[MAX_TAXA][MAX_TAXA], int num_taxa, double TD_arr[MAX_TAXA], int& index1, int& index2) {
    double min_distance = INT_MAX;
    for (int i = 0; i < num_taxa; i++) {
        if(arr[i][0]!=-1) {
            for (int j = i + 1; j < num_taxa; j++) {
                if(arr[j][0]!=-1){
                    double D_star = (num_taxa - 2) * arr[i][j] - TD_arr[i] - TD_arr[j];
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

//Pending
void updateDistanceMatrix(double arr[MAX_TAXA][MAX_TAXA], int num_taxa, int min_index, int max_index) {
    for (int k = 0; k < num_taxa; k++) {
        if (k != min_index && k != max_index) {
            arr[max_index][k] = ( arr[min_index][k] + arr[max_index][k] - arr[min_index][max_index]) / 2;
            arr[k][max_index] = arr[max_index][k];
        }
    }
    arr[min_index][0] = arr[0][min_index] = -1;
}

int main() {
    
    string file_name = "./examples/evolution.in";
    double arr[MAX_TAXA][MAX_TAXA];
    char seq[MAX_TAXA];
    Node* nodes[MAX_TAXA];
    //int num_taxa = read_DM_file(arr, seq, file_name, nodes);
    int num_taxa = readFromFile(arr, seq, file_name, nodes);
    printDistanceMatrix(arr, num_taxa, nodes);

    double TD_arr[MAX_TAXA];
    auto start_time = high_resolution_clock::now();

    for(int i=0 ; i<num_taxa -2; i++) {
        //{Need to plug in neighbor joining below to get indexes
        int index1;
        int index2;
        int n = num_taxa - i;
        //getRandomIndexes(index1, index2, n);
        totalDistance(arr, num_taxa, TD_arr);
        printTDMatrix(TD_arr, num_taxa);
        find_closest_pair(arr,num_taxa, TD_arr, index1, index2);
        //}
        
        int min_index = min(index1, index2);
        int max_index = max(index1, index2);
        double delta_ij = (TD_arr[min_index] - TD_arr[max_index]) / (n-2);
        double limb_length_i = (arr[min_index][max_index] + delta_ij)/2.0;
        double limb_length_j = (arr[min_index][max_index] - delta_ij)/2.0;
        updateDistanceMatrix(arr,num_taxa, min_index, max_index);
        string new_node_name = "(" + nodes[min_index]->node_name + nodes[max_index]->node_name + ")";
        cout<<new_node_name<<endl;
        Node* temp = new Node(new_node_name, nodes[min_index], nodes[max_index], limb_length_i, limb_length_j );
        nodes[max_index] = temp;
        nodes[min_index] = nullptr;
        //printDistanceMatrix(arr, num_taxa, nodes);
    }
    
    auto end_time = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end_time - start_time);
    printf("### \n Elapsed Time %" PRId64 "\n###\n", duration.count());


    int final_index1 = -1;
    int final_index2 = -1;

    for(int i=0 ; i<num_taxa ; i++) {
        if(arr[i][0]!=-1)
        {
            if(final_index1==-1)
                final_index1 = i;
            else
                final_index2 = i;
        }
    } 

    string root_node_name = "(" + nodes[final_index1]->node_name + nodes[final_index2]->node_name + ")";
    cout<<root_node_name<<endl;
    Node* root = new Node(root_node_name, nodes[final_index1], nodes[final_index2], arr[final_index1][final_index2]/2.0, arr[final_index1][final_index2]/2.0 );

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
