import numpy as np

# Load the .npz file
data = np.load('results/PS_impedance/transmission.npz')

# Print the names of the arrays in the file
print("Arrays in the file:", data.files)

# Iterate over each array and print its content
for array_name in data.files:
    print(f"\nContent of array '{array_name}':")
    print(data[array_name])