import pandas as pd

class BaseFormatter:
    def __init__(self, row):
        self.row = row

    @property
    def atom_type(self):
        return "{:4s}".format(self.row.get('atom_type', ''))

    @property
    def x(self):
        return "{:8.3f}".format(self.row.get('x', 0.0))

    @property
    def y(self):
        return "{:8.3f}".format(self.row.get('y', 0.0))

    @property
    def z(self):
        return "{:8.3f}".format(self.row.get('z', 0.0))

    @property
    def atom_name(self):
        return "{:4s}".format(self.row.get('atom_name', 1))
    
    @property
    def mol_name(self):
        return "{:3s}".format(self.row.get('mol_name', 'UNK'))

    @property
    def mol_num(self):
        return "{:4d}".format(self.row.get('mol_num', 1))
    
    def format_row(self):
        raise NotImplementedError("This method should be implemented in a child class")

class PDBFormatter(BaseFormatter):
    @property
    def atom_id(self):
        return "{:5d}".format(self.row.get('atom_id', 1))

    @property
    def atom_name(self):
        return "{:<4}".format(self.row.get('atom_name', '    '))

    @property
    def alt_loc(self):
        return "{:1}".format(self.row.get('alt_loc', ' '))

    @property
    def res_name(self):
        return "{:>3}".format(self.row.get('res_name', '   '))

    @property
    def chain_id(self):
        return "{:1}".format(self.row.get('chain_id', 'A'))

    @property
    def res_seq_num(self):
        return "{:4d}".format(self.row.get('res_seq_num', 1))

    @property
    def insertion_code(self):
        return "{:1}".format(self.row.get('insertion_code', ' '))

    @property
    def segment_id(self):
        return "{:<4}".format(self.row.get('segment_id', '    '))

    @property
    def element_symbol(self):
        return "{:>2}".format(self.row.get('element_symbol', '  '))

    @property
    def charge(self):
        return "{:<2}".format(self.row.get('charge', '  '))

    @property
    def occupancy(self):
        return "{:6.2f}".format(self.row.get('occupancy', 1.00))

    @property
    def temp_factor(self):
        return "{:6.2f}".format(self.row.get('temp_factor', 0.00))

    def format_row(self):
        return "ATOM  {} {}{}{} {}{}{}   {}{}{}  {}{}          {}{}".format(
            self.atom_id, self.atom_name, self.alt_loc, self.res_name, 
            self.chain_id, self.res_seq_num, self.insertion_code, 
            self.x, self.y, self.z, self.occupancy, self.temp_factor, 
            self.segment_id, self.element_symbol, self.charge)



class XYZFormatter(BaseFormatter):
    def format_row(self):
        return "{} {} {} {}".format(
            self.atom_name, self.x, self.y, self.z)


class DataFrameConverter:
    def __init__(self, formatter_class):
        self.formatter_class = formatter_class

    def convert(self, df):
        formatted_rows = [self.formatter_class(row).format_row() for _, row in df.iterrows()]
        return '\n'.join(formatted_rows)