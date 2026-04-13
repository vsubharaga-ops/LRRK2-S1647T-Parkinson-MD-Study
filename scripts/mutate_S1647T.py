from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa

class MutantSelect(Select):
    def accept_residue(self, residue):
        return is_aa(residue)

parser = PDBParser(QUIET=True)
structure = parser.get_structure("lrrk2", "your_lrrk2.pdb")

mutated = False
for model in structure:
    for chain in model:
        for residue in chain:
            if residue.get_id()[1] == 1647:
                if residue.get_resname() == "SER":
                    residue.resname = "THR"
                    print(f"Mutated SER 1647 → THR in chain {chain.id}")
                    mutated = True

if not mutated:
    print("WARNING: Residue 1647 not found or not SER — check your PDB numbering!")

io = PDBIO()
io.set_structure(structure)
io.save("mutant.pdb", MutantSelect())
print("Saved mutant.pdb")
