# In-Place Loop Subdivision using Half-Edge Data Structure

This code uses the Half-Edge Data Structure [Bischoff et. al (2002)] as an efficient mesh traversal method to perform in-place loop subdivision [Loop (1987)] to minimize memory utilization for the subdivision process. The algorithm uses edge split and edge flip operations to generate subdivision vertices without operating on a separate mesh.

Bischoff, B. S., Botsch, M., Steinberg, S., Bischoff, S., Kobbelt, L., & Aachen, R. (2002). OpenMesh–a generic and efficient polygon mesh data structure. In In openSG symposium (Vol. 18).
Loop, C. (1987). Smooth subdivision surfaces based on triangles.

## Setup and Usage

Run the following commands in your bash terminal

1. cd ./Code/Subdivision/
2. make
3. Subdivision <input_file_path> <output_file_path>

Both input and output files should be .obj file formats.
