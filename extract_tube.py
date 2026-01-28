"""
Extract height, outer diameter, and wall thickness from a hollow cylinder (tube) STEP file.
"""

from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.GeomAbs import GeomAbs_Cylinder
from OCC.Core.BRepBndLib import brepbndlib
from OCC.Core.Bnd import Bnd_Box


def extract_tube_dimensions(filepath):
    """Extract height, outer diameter, and thickness from a tube STEP file."""
    
    # Load STEP file
    reader = STEPControl_Reader()
    if reader.ReadFile(filepath) != 1:
        raise Exception(f"Error reading STEP file: {filepath}")
    reader.TransferRoots()
    shape = reader.OneShape()
    
    # Find all cylindrical surfaces and get their radii
    radii = []
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    while explorer.More():
        face = explorer.Current()
        adaptor = BRepAdaptor_Surface(face)
        if adaptor.GetType() == GeomAbs_Cylinder:
            radii.append(adaptor.Cylinder().Radius())
        explorer.Next()
    
    if len(radii) < 2:
        raise Exception("Not a hollow cylinder - need at least 2 cylindrical surfaces")
    
    # Get unique radii (inner and outer)
    unique_radii = sorted(set(round(r, 6) for r in radii))
    inner_radius = min(unique_radii)
    outer_radius = max(unique_radii)
    
    # Get height from bounding box
    bbox = Bnd_Box()
    brepbndlib.Add(shape, bbox)
    xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
    height = zmax - zmin
    
    return {
        "height": round(height, 4),
        "outer_diameter": round(outer_radius * 2, 4),
        "inner_diameter": round(inner_radius * 2, 4),
        "thickness": round(outer_radius - inner_radius, 4)
    }


MODEL_PATH = "model.STEP"


if __name__ == "__main__":
    dims = extract_tube_dimensions(MODEL_PATH)
    
    print(f"Height:         {dims['height']} mm")
    print(f"Outer Diameter: {dims['outer_diameter']} mm")
    print(f"Inner Diameter: {dims['inner_diameter']} mm")
    print(f"Wall Thickness: {dims['thickness']} mm")
