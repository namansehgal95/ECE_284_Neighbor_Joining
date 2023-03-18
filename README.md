# ECE_284_Neighbor_Joining

Run the code using the command:

    g++ neighbor_join.cpp -o a.exe; ./a.exe

Following this, run Graphviz using the below command (g.gv dumped from the above command)

    dot -Tpng g.gv -o file.png

This will generate a png file called file.png which can be viewed with any image viewer. 

Note: Will need to download Graphviz. It can be downloaded from https://graphviz.org/download/

For running the GPU code, please compile using the following command:

    nvcc gpu_neighbor_joining.cu -o gpu_nj
    
Then the executable can be run by:

    ./gpu_nj
    
