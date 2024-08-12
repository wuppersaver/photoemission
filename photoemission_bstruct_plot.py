from pathlib import Path
from optados_output_class import OptaDOSOutput

class BandStructPlot2D:

    def __init__(self) -> None:
        pass
    

    def from_files(self,path:Path):

        return self

    def set_contributions(self,path:Path):
        pass
    
    def set_structure(self,path:Path):
        
        pass

    def set_bands(self,path:Path):
        if self.structure == None:
            self.set_structure(path.parent)
        
        pass
    