import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Function to read data from the C++ generated text file
def read_data(file_path):
    with open(file_path, 'r') as file:
        # Skip the first line
        next(file)
        data = []
        for line in file:
            try:
                parts = line.split()
                x_index = parts.index('x:')
                y_index = parts.index('y')
                z_index = parts.index('z')
                row = [float(parts[x_index + 1]), float(parts[y_index + 2]), float(parts[z_index + 2])]
                data.append(row)
            except (ValueError, IndexError):
                print(f"Skipping line: {line.strip()} (non-numeric data or unexpected format)")

    return list(zip(*data))


# Path to the C++ generated output file
file_path = r'D:\ODE Solver\ODE solver\output.txt'

# Read data from the file
x, y, z = read_data(file_path)

# Plot the results in 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x, y, z, label='Lorenz System')
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Z-axis')
ax.set_title('Lorenz System')
ax.legend()

# Show the plot
plt.show()
