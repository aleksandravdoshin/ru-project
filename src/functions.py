import pandas as pd 
import numpy as np
from scipy.spatial.transform import Rotation as R

def read_pdb(filename):
    atoms = []
    coords = []
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atoms.append(line[76:78].strip())
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])

    df = pd.DataFrame(coords, columns=["x", "y", "z"])
    df['atom'] = atoms
    return df


def read_cif_data(filename):
    parameters = {}
    atom_data = []
    
    # Флаг, указывающий на то, что мы начали читать данные атомов
    reading_atoms = False
    
    with open(filename, 'r') as file:
        lines = file.readlines()
        for line in lines:
            line = line.strip()
            if "_cell_length_a" in line:
                parameters["a"] = float(line.split()[1])
            elif "_cell_length_b" in line:
                parameters["b"] = float(line.split()[1])
            elif "_cell_length_c" in line:
                parameters["c"] = float(line.split()[1])
            elif "_cell_angle_alpha" in line:
                parameters["alpha"] = float(line.split()[1])
            elif "_cell_angle_beta" in line:
                parameters["beta"] = float(line.split()[1])
            elif "_cell_angle_gamma" in line:
                parameters["gamma"] = float(line.split()[1])
            elif "_atom_site_type_symbol" in line:
                reading_atoms = True
                continue
            elif reading_atoms and line.startswith("loop_"):
                reading_atoms = False
            elif reading_atoms and len(line.split()) == 7:
                parts = line.split()
                atom_data.append([parts[0], float(parts[3]), float(parts[4]), float(parts[5])])
    
    df = pd.DataFrame(atom_data, columns=["atom", "x", "y", "z"])
    
    return parameters, df


def write_cif(filename, parameters, df):
    with open(filename, 'w') as file:
        file.write("# generated using our script\n")
        file.write("data_structure\n")
        file.write("_symmetry_space_group_name_H-M   'P 1'\n")
        file.write(f"_cell_length_a   {parameters['a']}\n")
        file.write(f"_cell_length_b   {parameters['b']}\n")
        file.write(f"_cell_length_c   {parameters['c']}\n")
        file.write(f"_cell_angle_alpha   {parameters['alpha']}\n")
        file.write(f"_cell_angle_beta   {parameters['beta']}\n")
        file.write(f"_cell_angle_gamma   {parameters['gamma']}\n")
        file.write("_symmetry_Int_Tables_number   1\n")
        file.write("_chemical_formula_structural   formula_here\n")  # Замените formula_here на вашу формулу, если она у вас есть
        file.write("_chemical_formula_sum   'formula_sum_here'\n")  # Замените formula_sum_here на суммарную формулу, если она у вас есть
        file.write("_cell_volume   volume_here\n")  # Замените volume_here на объем ячейки, если он у вас есть
        file.write("_cell_formula_units_Z   units_here\n")  # Замените units_here на количество формульных единиц в ячейке, если это известно
        file.write("loop_\n")
        file.write(" _symmetry_equiv_pos_site_id\n")
        file.write(" _symmetry_equiv_pos_as_xyz\n")
        file.write("  1  'x, y, z'\n")
        file.write("loop_\n")
        file.write(" _atom_site_type_symbol\n")
        file.write(" _atom_site_label\n")
        file.write(" _atom_site_symmetry_multiplicity\n")
        file.write(" _atom_site_fract_x\n")
        file.write(" _atom_site_fract_y\n")
        file.write(" _atom_site_fract_z\n")
        file.write(" _atom_site_occupancy\n")

        for index, row in df.iterrows():
            file.write(f"  {row['atom']}  {row['atom']}{index+1}  1  {row['x']}  {row['y']}  {row['z']}  1.0\n")


def lattice_to_cartesian(coord, lattice_params):
    """
    Converts fractional coordinates to Cartesian coordinates using lattice parameters.
    Corrects the transformation matrix for non-orthogonal systems.
    """
    a = lattice_params['a']
    b = lattice_params['b']
    c = lattice_params['c']
    alpha = np.deg2rad(lattice_params['alpha'])
    beta = np.deg2rad(lattice_params['beta'])
    gamma = np.deg2rad(lattice_params['gamma'])
    
    # Vectors of the lattice in Cartesian coordinates
    ax = a
    ay = 0
    az = 0

    bx = b * np.cos(gamma)
    by = b * np.sin(gamma)
    bz = 0

    cx = c * np.cos(beta)
    cy = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    cz = c * np.sqrt(1 - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2 + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma)) / np.sin(gamma)

    # Transformation matrix
    transform_matrix = np.array([
        [ax, bx, cx],
        [ay, by, cy],
        [az, bz, cz]
    ])

    return np.dot(coord, transform_matrix.T)


def write_xyz(df, filename="output.xyz"):
    with open(filename, "w") as file:
        file.write(f"{len(df)}\n")
        file.write("Generated using our script\n")
        for index, row in df.iterrows():
            file.write(f"{row['atom']} {row['x']} {row['y']} {row['z']}\n")


def cartesian_to_lattice(cartesian_coords, lattice_params):
    """
    Converts Cartesian coordinates to fractional coordinates using lattice parameters.
    Uses the inverse of the transformation matrix used for the original conversion.
    """
    a = lattice_params['a']
    b = lattice_params['b']
    c = lattice_params['c']
    alpha = np.deg2rad(lattice_params['alpha'])
    beta = np.deg2rad(lattice_params['beta'])
    gamma = np.deg2rad(lattice_params['gamma'])

    # Vectors of the lattice in Cartesian coordinates
    ax = a
    ay = 0
    az = 0

    bx = b * np.cos(gamma)
    by = b * np.sin(gamma)
    bz = 0

    cx = c * np.cos(beta)
    cy = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    cz = c * np.sqrt(1 - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2 + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma)) / np.sin(gamma)

    # Transformation matrix
    transform_matrix = np.array([
        [ax, bx, cx],
        [ay, by, cy],
        [az, bz, cz]
    ])

    # Inverse of the transformation matrix
    inverse_transform_matrix = np.linalg.inv(transform_matrix)

    # Convert to fractional coordinates
    return np.dot(cartesian_coords, inverse_transform_matrix.T)


def rotate_dataframe_euler(df):
    """
    Rotate a DataFrame with xyz coordinates using Euler angles.

    Parameters:
    df (pandas.DataFrame): DataFrame containing the coordinates of points in columns 'x', 'y', 'z'.

    Returns:
    pandas.DataFrame: The rotated DataFrame.
    """
    # Calculate the geometric center
    gc_initial = df[['x', 'y', 'z']].mean()

    # Move GC to origin
    df_centered = df[['x', 'y', 'z']] - gc_initial

    # Generate random Euler angles (in radians)
    angles = np.random.rand(3) * 2 * np.pi  # Random angles for rotation around x, y, z axes

    # Create a rotation object based on Euler angles
    rotation = R.from_euler('xyz', angles)

    # Apply the rotation
    rotated_points = rotation.apply(df_centered)

    # Move GC back to initial position
    df_rotated = pd.DataFrame(rotated_points, columns=['x', 'y', 'z']) + gc_initial

    return df_rotated

