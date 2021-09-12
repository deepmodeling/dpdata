import dpdata.amber.md
import dpdata.amber.sqm
from dpdata.format import Format


@Format.register("amber/md")
class AmberMDFormat(Format):
    def from_system(self, file_name=None, parm7_file=None, nc_file=None, use_element_symbols=None, **kwargs):
        # assume the prefix is the same if the spefic name is not given
        if parm7_file is None:
            parm7_file = file_name + ".parm7"
        if nc_file is None:
            nc_file = file_name + ".nc"
        return dpdata.amber.md.read_amber_traj(parm7_file=parm7_file, nc_file=nc_file, use_element_symbols=use_element_symbols, labeled=False)

    def from_labeled_system(self, file_name=None, parm7_file=None, nc_file=None, mdfrc_file=None, mden_file=None, mdout_file=None, use_element_symbols=None, **kwargs):
        # assume the prefix is the same if the spefic name is not given
        if parm7_file is None:
            parm7_file = file_name + ".parm7"
        if nc_file is None:
            nc_file = file_name + ".nc"
        if mdfrc_file is None:
            mdfrc_file = file_name + ".mdfrc"
        if mden_file is None:
            mden_file = file_name + ".mden"
        if mdout_file is None:
            mdout_file = file_name + ".mdout"
        return dpdata.amber.md.read_amber_traj(parm7_file, nc_file, mdfrc_file, mden_file, mdout_file, use_element_symbols)


@Format.register("sqm/out")
class SQMOutFormat(Format):
    def from_system(self, fname, **kwargs):
        '''
        Read from ambertools sqm.out
        '''
        return dpdata.amber.sqm.parse_sqm_out(fname)
    
    def from_labeled_system(self, fname, **kwargs):
        '''
        Read from ambertools sqm.out
        '''
        data = dpdata.amber.sqm.parse_sqm_out(fname)
        assert "forces" in list(data.keys()), f"No forces in {fname}"
        return data

@Format.register("sqm/in")
class SQMINFormat(Format):
    def to_system(self, data, fname=None, frame_idx=0, **kwargs):
        """
        Generate input files for semi-emperical calculation in sqm software
        """
        return dpdata.amber.sqm.make_sqm_in(data, fname, frame_idx, **kwargs)
