import pandas as pd
import re
import duckdb

def replace_header_components(header_str, atoms, bonds, angles, dihedrals, impropers):
    """
    Replace specific lines in the given header string with provided values.
    
    :param header_str: The input string
    :param atoms: New value for atoms
    :param bonds: New value for bonds
    :param angles: New value for angles
    :param dihedrals: New value for dihedrals
    :param impropers: New value for impropers
    :return: Modified header string
    """
    
    replacements = {
        "atoms": atoms,
        "bonds": bonds,
        "angles": angles,
        "dihedrals": dihedrals,
        "impropers": impropers
    }
    
    for key, value in replacements.items():
        if value is None:
            continue
        header_str = re.sub(f"(\\d+)\\s+{key}", f"{value:>12} {key}", header_str)
    header_str = '\n'.join([i.strip() for i in header_str.split('\n')])
    return header_str

class LammpsData:
    def __init__(self, file_path):
        with open(file_path, 'r') as f:
            self.lines = f.readlines()
        self.header = None
        self.atoms = None
        self.bonds = None
        self.angles = None
        self.dihedrals = None
        self.impropers = None

        self.pair_coeffs = None
        self.masses_coeffs = None
        self.bonds_coeffs = None
        self.angles_coeffs = None
        self.dihedrals_coeffs = None
        self.impropers_coeffs = None


    def _find_section_lines(self, section):
        start_key = section.strip()
        try:
            start_index = [i.strip()[:len(start_key)]
                           for i in self.lines].index(start_key) + 2
        except ValueError:
            return None, None

        end_index = start_index
        while end_index < len(self.lines) and self.lines[end_index].strip():
            end_index += 1

        return start_index, end_index


    def _create_dataframe_from_section(self, section, columns, dtypes=None):
        start_index, end_index = self._find_section_lines(section)
        if start_index is None:
            return pd.DataFrame(columns=columns)

        data = []
        for line in self.lines[start_index:end_index]:
            parts = line.split('#')
            values = parts[0].split()
            if len(parts) > 1:
                comment = parts[1].strip()
                values.append(comment)
            data.append(values)
        
        df = pd.DataFrame(data, columns=columns)
        if dtypes:
            df = df.astype(dtypes)
        return df


    def _process_coefficient_section(self, section):
        start_index, end_index = self._find_section_lines(section)
        if start_index is None:
            return {}

        coeff_data = {}
        for line in self.lines[start_index:end_index]:
            # Splitting line into components and comment
            components, comment = line.split("#")

            # Extracting ID, style, and coefficients
            components_list = components.split()
            coeff_id = int(components_list[0])
            style = components_list[1]
            coefficients = [float(val) for val in components_list[2:]]

            # Storing in the dictionary
            coeff_data[coeff_id] = {
                'style': style,
                'coefficients': coefficients,
                'comment': comment.strip()
            }

        return coeff_data


    def _process_header(self):
        # Check where the header ends; typically before "Atoms" section or another main section.
        end_index = self._find_section_lines('Atoms')[0] - 2
        header_lines = self.lines[:end_index]

        return ''.join(header_lines)


    def _process_masses(self, section):
        start_index, end_index = self._find_section_lines(section)
        if start_index is None:
            return None
    

    def read(self):
        # Atoms
        atom_dtypes = {
            'id': int,
            'mol_id': int,
            'atom_type': int,
            'q': float,
            'x': float,
            'y': float,
            'z': float
        }
        try:
            self.masses = self._create_dataframe_from_section('Masses',
                                                            ['id', 'mass', 
                                                                'atom'],
                                                            {'id': int, 'mass': float})
        except Exception:
            self.masses = None
            print('No masses section found')

        try:
            self.pair_coeffs = self._create_dataframe_from_section('Pair Coeffs',
                                                            ['id', 'epselon', 
                                                              'sigma',  'atom'],
                                                              {
                                                                'id': int,
                                                                'epselon': float,
                                                                'sigma': float
                                                            })
            
        except Exception:
            self.pair_coeffs = None
            print('No pair_coeffs section found')
        
        try:
            self.bonds_coeffs = self._create_dataframe_from_section('Bond Coeffs',
                                                                    ['id', 'stiffness', 'r', 'comment'],
                                                                    {'id': int,
                                                                        'stiffness': float,
                                                                        'r': float})
            self.bonds_coeffs['atom0'] = self.bonds_coeffs.comment.apply(lambda x: x.split()[0])
            self.bonds_coeffs['atom1'] = self.bonds_coeffs.comment.apply(lambda x: x.split()[1])
            
            self.bonds_coeffs = self.bonds_coeffs.drop(columns={'comment'})

        except Exception:
            self.bonds_coeffs = None
            print('No bonds_coeffs section found')

        try:
            self.angles_coeffs = self._process_coefficient_section('Angle Coeffs') 
        except Exception:
            self.angles_coeffs = None
            print('No angles_coeffs section found')
        
        try:
            self.dihedrals_coeffs = self._process_coefficient_section('Dihedral Coeffs')
        except Exception:
            self.dihedrals_coeffs = None
            print('No dihedrals_coeffs section found')

        try:
            self.impropers_coeffs = self._process_coefficient_section('Improper Coeffs')
        except Exception:
            self.impropers_coeffs = None
            print('No impropers_coeffs section found')    
        
        try:
            self.atoms = self._create_dataframe_from_section('Atoms',
                                                            ['id', 'mol_id', 'atom_type',
                                                                'q', 'x', 'y', 'z'],
                                                            atom_dtypes)
        except Exception:
            self.atoms = None
            print('No atoms section found')        

        try:
            # Bonds
            self.bonds = self._create_dataframe_from_section('Bonds',
                                                            ['id', 'bond_type',
                                                                'atom1_id', 'atom2_id'],
                                                            int)
        except Exception:
            self.bonds = None
            print('No bonds section found')
        
        try:
            # Angles
            self.angles = self._create_dataframe_from_section('Angles',
                                                            ['id', 'angle_type', 'atom1_id',
                                                                'atom2_id', 'atom3_id'],
                                                            int)
        except Exception:
            self.angles = None
            print('No angles section found')

        
        try:
            # Dihedrals
            self.dihedrals = self._create_dataframe_from_section('Dihedrals',
                                                                ['id', 'dihedral_type', 'atom1_id',
                                                                    'atom2_id', 'atom3_id', 'atom4_id'],
                                                                int)
        except Exception:
            self.dihedrals = None
            print('No dihedrals section found')

        try:
            # # Impropers
            self.impropers = self._create_dataframe_from_section('Impropers',
                                                                ['id', 'improper_type', 'atom1_id',
                                                                    'atom2_id', 'atom3_id', 'atom4_id'],
                                                                int)
        except Exception:
            self.impropers = None
            print('No impropers section found')


        try:
            self.header = self._process_header()
        except Exception:
            self.header = None
            print('No header found')


    def choose_mol(self, id):
        self.mol_atoms = self.atoms[self.atoms.mol_id == id].copy()
        alken_ids = set(
            self.mol_atoms.id.values)

        self.mol_bonds = self.bonds[self.bonds.atom1_id.isin(alken_ids) &
                        self.bonds.atom2_id.isin(alken_ids)].copy()

        self.mol_angles = self.angles[self.angles.atom1_id.isin(alken_ids) &
                        self.angles.atom2_id.isin(alken_ids) &
                        self.angles.atom3_id.isin(alken_ids)].copy()

        self.mol_dihedrals = self.dihedrals[self.dihedrals.atom1_id.isin(alken_ids) &
                            self.dihedrals.atom2_id.isin(alken_ids) &
                            self.dihedrals.atom3_id.isin(alken_ids) &
                            self.dihedrals.atom4_id.isin(alken_ids)].copy()

        self.mol_impropers = self.impropers[self.impropers.atom1_id.isin(alken_ids) &
                            self.impropers.atom2_id.isin(alken_ids) &
                            self.impropers.atom3_id.isin(alken_ids) &
                            self.impropers.atom4_id.isin(alken_ids)].copy()


    def append_mol(self):
        delta_id = self.atoms['id'].max() - \
            self.mol_atoms['id'].min() + 1
        
        self.mol_atoms['id'] += delta_id
        self.mol_bonds['atom1_id'] += delta_id
        self.mol_bonds['atom2_id'] += delta_id
        if self.mol_angles is not None:
            self.mol_angles['atom1_id'] += delta_id
            self.mol_angles['atom2_id'] += delta_id
            self.mol_angles['atom3_id'] += delta_id
        if self.mol_dihedrals is not None:
            self.mol_dihedrals['atom1_id'] += delta_id
            self.mol_dihedrals['atom2_id'] += delta_id
            self.mol_dihedrals['atom3_id'] += delta_id
            self.mol_dihedrals['atom4_id'] += delta_id
        if self.mol_impropers is not None:
            self.mol_impropers['atom1_id'] += delta_id
            self.mol_impropers['atom2_id'] += delta_id
            self.mol_impropers['atom3_id'] += delta_id
            self.mol_impropers['atom4_id'] += delta_id
        
        self.mol_atoms['mol_id'] = self.atoms.mol_id.max() + 1
        

        self.mol_bonds['id'] += 1 + self.bonds['id'].max() - self.mol_bonds['id'].min()
        if self.mol_angles is not None:
            self.mol_angles['id'] += 1 + self.angles['id'].max() - self.mol_angles['id'].min()
        if self.mol_dihedrals is not None:
            self.mol_dihedrals['id'] += 1 + self.dihedrals['id'].max() - self.mol_dihedrals['id'].min()
        if self.mol_impropers is not None:
            self.mol_impropers['id'] += 1 + self.impropers['id'].max() - self.mol_impropers['id'].min()
        
        self.atoms = pd.concat([self.atoms, self.mol_atoms]).reset_index(drop=True)
        self.bonds = pd.concat([self.bonds, self.mol_bonds]).reset_index(drop=True)
        if self.mol_angles is not None:
            self.angles = pd.concat([self.angles, self.mol_angles]).reset_index(drop=True)
        if self.mol_dihedrals is not None:
            self.dihedrals = pd.concat([self.dihedrals, self.mol_dihedrals]).reset_index(drop=True)
        if self.impropers is not None:
            self.impropers = pd.concat([self.impropers, self.mol_impropers]).reset_index(drop=True)

        len_angles = None if self.angles is None else len(self.angles)
        len_dihedrals = None if self.dihedrals is None else len(self.dihedrals)
        len_imp = None if self.impropers is None else len(self.impropers)


        self.header = replace_header_components(self.header, len(self.atoms), 
                                                len(self.bonds), len_angles, 
                                                len_dihedrals, len_imp)
        

    
    @property
    def atoms_map(self):
        return {i: j['comment'] for i, j in self.masses.items()}

    @property
    def xyzd(self):
        
        xyzd = self.atoms.copy()
        xyzd['atom_type'] = xyzd['atom_type'].map(self.atoms_map)
        return xyzd
    
    @property
    def bonds_distances(self):
        connection = duckdb.connect(database=':memory:', read_only=False)
        connection.register('df', self.atoms)
        connection.register('bonds', self.bonds)
        s = '''select bonds.*, 
                            ((df1.x - df2.x) ** 2 + (df1.y - df2.y) ** 2 + (df1.z - df2.z) ** 2) ** 0.5 as distance
                                from bonds
                left join df df1 on df1.id = bonds.atom1_id
                left join df df2 on df2.id = bonds.atom2_id'''
        
        bonds = connection.execute(s).fetchdf()
        return bonds
    
    @property
    def num_neigbours(self):
        connection = duckdb.connect(database=':memory:', read_only=False)
        connection.register('dft', self.atoms)
        connection.register('bonds', self.bonds)
        connection.register('xyz', self.xyzd)
        query = """
        WITH BondCounts AS (
            SELECT atom1_id AS atom_id FROM bonds
            UNION ALL
            SELECT atom2_id AS atom_id FROM bonds
        )
        SELECT dft.id, dft.atom_type, COUNT(BondCounts.atom_id) AS bond_count
        FROM dft
        LEFT JOIN BondCounts ON dft.id = BondCounts.atom_id
        GROUP BY dft.id, dft.atom_type
        ORDER BY dft.id
        """
        df = connection.execute(query).fetchdf()
        query = """
        select df.id, xyz.atom_type, mol_id, q, bond_count from df
        left join xyz on df.id = xyz.id
        """
        df = connection.execute(query).fetchdf()
        return df 
    
    @property
    def list_of_neighbours(self):
        neighbours = self.bonds_distances.copy()[['atom1_id', 'atom2_id', 'distance', 'bond_type']]
        atom_map = self.xyzd[['id', 'atom_type']].set_index('id')['atom_type'].to_dict()
        neighbours['atom1'] = neighbours.atom1_id.map(atom_map)
        neighbours['atom2'] = neighbours.atom2_id.map(atom_map)
        df = pd.concat([neighbours[['atom1_id', 'atom2_id', 'distance', 'atom1', 'atom2', 'bond_type']],
                        neighbours[['atom2_id', 'atom1_id', 'distance', 'atom2', 'atom1', 'bond_type']]\
                            .rename(columns={'atom2_id': 'atom1_id', 'atom1_id': 'atom2_id', 
                                             'atom2': 'atom1', 'atom1': 'atom2'})])
        return df

        
class LammpsDataWriter:

    def __init__(self, reader_instance):
        self.data = reader_instance

    def write_header(self, file):
        file.write(f"{self.data.header}\n\n")

    def write_masses(self, file):
        file.write("Masses\n\n")
        self.data.masses = self.data.masses.astype({'id': 'int', 'mass': 'float'})
        for _, row in self.data.masses.iterrows():
            file.write(f"   {row['type']} {row['mass']} # {row['comment']}\n")
        file.write("\n")

    def write_coeffs(self, file, section_name, data_attr):
        file.write(f"{section_name}\n\n")
        coeffs_data = getattr(self.data, data_attr)
        for key, val in coeffs_data.items():
            coeffs = ' '.join(map(str, val['coefficients']))
            file.write(
                f"   {key} {val['style']} {coeffs} # {val['comment']}\n")
        file.write("\n")

    def write_atoms(self, file):
        file.write("Atoms\n\n")
        self.data.atoms = self.data.atoms.astype({'id': 'int', 'mol_id': 'int', 'atom_type': 'int', 'q': 'float'})
        for _, row in self.data.atoms.iterrows():
            file.write(
                f"{int(row['id']):>8} {int(row['mol_id']):>8} {int(row['atom_type']):>8} {row['q']:>12.5f} {row['x']:>12.5f} {row['y']:>12.5f} {row['z']:>12.5f}\n")
        file.write("\n")

    def write_bonds(self, file):
        file.write("Bonds\n\n")
        self.data.bonds = self.data.bonds.astype({'id': 'int', 'bond_type': 'int', 'atom1_id': 'int', 'atom2_id': 'int'})
        for _, row in self.data.bonds.iterrows():
            file.write(
                f"{row['id']:>8} {row['bond_type']:>8} {row['atom1_id']:>8} {row['atom2_id']:>8}\n")
        file.write("\n")

    def write_angles(self, file):
        file.write("Angles\n\n")
        self.data.angles = self.data.angles.astype({'id': 'int', 'angle_type': 'int', 'atom1_id': 'int', 'atom2_id': 'int', 'atom3_id': 'int'})
        for _, row in self.data.angles.iterrows():
            file.write(
                f"{row['id']:>8} {row['angle_type']:>8} {row['atom1_id']:>8} {row['atom2_id']:>8} {row['atom3_id']:>8}\n")
        file.write("\n")

    def write_dihedrals(self, file):
        file.write("Dihedrals\n\n")
        self.data.dihedrals = \
        self.data.dihedrals.astype({'id': 'int', 'dihedral_type': 'int', 
        'atom1_id': 'int', 'atom2_id': 'int', 
        'atom3_id': 'int', 'atom4_id': 'int'})
        for _, row in self.data.dihedrals.iterrows():
            file.write(
                f"{row['id']:>8} {row['dihedral_type']:>8} {row['atom1_id']:>8} {row['atom2_id']:>8} {row['atom3_id']:>8} {row['atom4_id']:>8}\n")
        file.write("\n")

    def write_impropers(self, file):
        file.write("Impropers\n\n")
        self.data.impropers = \
        self.data.impropers.astype({'id': 'int', 'improper_type': 'int', 
        'atom1_id': 'int', 'atom2_id': 'int', 
        'atom3_id': 'int', 'atom4_id': 'int'})
        if self.data.impropers is None:
            return
        for _, row in self.data.impropers.iterrows():
            file.write(
                f"{row['id']:>8} {row['improper_type']:>8} {row['atom1_id']:>8} {row['atom2_id']:>8} {row['atom3_id']:>8} {row['atom4_id']:>8}\n")
        file.write("\n")

    def write_to_file(self, file_path):
        with open(file_path, 'w') as f:
            self.write_header(f)
            self.write_atoms(f)
            self.write_bonds(f)
            if self.data.angles is not None:
                self.write_angles(f)
            if self.data.dihedrals is not None:
                self.write_dihedrals(f)
            if self.data.impropers is not None:
                self.write_impropers(f)
