# ECE_284_Neighbor_Joining

Run the code using the command:

    g++ neighbor_join.cpp -o cpu_nj
       
The code can be executed by the command:
    
    ./cpu_nj

For running the GPU code, please compile using the following command:

    nvcc gpu_neighbor_joining.cu -o gpu_nj
    
Then the executable can be run by:

    ./gpu_nj
    

For running the gpu code you only need a gpu that runs successfully. We used dsmlp platform and the GPU availalbe to us there. 
We used a launch.sh script to create a pod with GPU instance and ran the results in our environment. We used the script for the course CSE 160 (as we couldn't launch a pod and hold on to it with the ECE 284 pod)
Probably the TA can also try to use our method which is:

    launch.sh -g 1 -s -i ucsdets/nvcr-cuda:latest -W CSE160_WI23_A00
    
Once the pod is launched the above command can be used to run the GPU code.
The input file can be modified on line #470 in the file gpu_neighbor_joining.cu. The current file that is selected is IN1000.in but any of the INxxx.in file can be selected. The final matrices of CPU and GPU were compared. The CPU input file can also be changed similarly on line #197.

This way you can setup the environment and reproduce the results.
