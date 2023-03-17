import numpy as np

num_taxa = 200
assert num_taxa > 1
# Generate a random matrix with values between 0.1 and 5
random_matrix = np.random.uniform(low=0.1, high=5, size=(num_taxa, num_taxa))

# Create a symmetric matrix from the upper triangular part of the random matrix
symmetric_matrix = np.triu(random_matrix) + np.triu(random_matrix, 1).T
np.fill_diagonal(symmetric_matrix, 0)
symmetric_matrix = np.round(symmetric_matrix, decimals=2)

print(random_matrix)
print(symmetric_matrix)

file = open("IN" + str(num_taxa) + ".in", "w")
file.write(str(num_taxa) + "\n")
file.write("A0")
for i in range(1,num_taxa):
    my_str ="\n"
    my_str += ("A" + str(i)).ljust(6)
    for j in range(0,i):
        my_str = my_str + " "  + str(symmetric_matrix[i][j]).ljust(4)
    file.write(my_str)
file.close()

import os
print (os.getcwd())