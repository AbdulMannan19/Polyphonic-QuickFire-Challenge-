"""
STEP File Analyzer
Extracts geometric and topological information from STEP files.
"""

from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.BRepBndLib import brepbndlib
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepGProp import brepgprop
from OCC.Core.GProp import GProp_GProps
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import (
    TopAbs_SOLID, TopAbs_FACE, TopAbs_EDGE, 
    TopAbs_VERTEX, TopAbs_SHELL, TopAbs_WIRE
)
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface, BRepAdaptor_Curve
from OCC.Core.GeomAbs import (
    GeomAbs_Plane, GeomAbs_Cylinder, GeomAbs_Cone,
    GeomAbs_Sphere, GeomAbs_Torus, GeomAbs_BSplineSurface,
    GeomAbs_Line, GeomAbs_Circle, GeomAbs_Ellipse,
    GeomAbs_BSplineCurve
)
from OCC.Core.BRep import BRep_Tool
import json
import math


def load_step(filepath):
    """Load a STEP file and return the shape."""
    reader = STEPControl_Reader()
    status = reader.ReadFile(filepath)
    if status != 1:
        raise Exception(f"Error reading STEP file: {filepath}")
    reader.TransferRoots()
    return reader.OneShape()


def get_bounding_box(shape):
    """Extract bounding box dimensions."""
    bbox = Bnd_Box()
    brepbndlib.Add(shape, bbox)
    xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
    return {
        "min_point": {"x": xmin, "y": ymin, "z": zmin},
        "max_point": {"x": xmax, "y": ymax, "z": zmax},
        "dimensions": {
            "width": xmax - xmin,
            "depth": ymax - ymin,
            "height": zmax - zmin
        },
        "center": {
            "x": (xmin + xmax) / 2,
            "y": (ymin + ymax) / 2,
            "z": (zmin + zmax) / 2
        }
    }


def get_mass_properties(shape):
    """Extract volume, surface area, and center of mass."""
    props = GProp_GProps()
    
    # Volume properties
    brepgprop.VolumeProperties(shape, props)
    volume = props.Mass()
    center_of_mass = props.CentreOfMass()
    
    # Surface properties
    brepgprop.SurfaceProperties(shape, props)
    surface_area = props.Mass()
    
    return {
        "volume": volume,
        "surface_area": surface_area,
        "center_of_mass": {
            "x": center_of_mass.X(),
            "y": center_of_mass.Y(),
            "z": center_of_mass.Z()
        }
    }


def get_topology_counts(shape):
    """Count topological entities (solids, faces, edges, vertices)."""
    counts = {
        "solids": 0,
        "shells": 0,
        "faces": 0,
        "wires": 0,
        "edges": 0,
        "vertices": 0
    }
    
    for topo_type, name in [
        (TopAbs_SOLID, "solids"),
        (TopAbs_SHELL, "shells"),
        (TopAbs_FACE, "faces"),
        (TopAbs_WIRE, "wires"),
        (TopAbs_EDGE, "edges"),
        (TopAbs_VERTEX, "vertices")
    ]:
        explorer = TopExp_Explorer(shape, topo_type)
        while explorer.More():
            counts[name] += 1
            explorer.Next()
    
    return counts


def get_surface_type_name(surface_type):
    """Convert surface type enum to readable name."""
    type_map = {
        GeomAbs_Plane: "plane",
        GeomAbs_Cylinder: "cylinder",
        GeomAbs_Cone: "cone",
        GeomAbs_Sphere: "sphere",
        GeomAbs_Torus: "torus",
        GeomAbs_BSplineSurface: "bspline_surface"
    }
    return type_map.get(surface_type, "other")


def get_curve_type_name(curve_type):
    """Convert curve type enum to readable name."""
    type_map = {
        GeomAbs_Line: "line",
        GeomAbs_Circle: "circle",
        GeomAbs_Ellipse: "ellipse",
        GeomAbs_BSplineCurve: "bspline_curve"
    }
    return type_map.get(curve_type, "other")


def analyze_faces(shape):
    """Analyze all faces and extract their geometric properties."""
    faces = []
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    face_id = 0
    
    while explorer.More():
        face = explorer.Current()
        adaptor = BRepAdaptor_Surface(face)
        surface_type = adaptor.GetType()
        
        face_info = {
            "id": face_id,
            "type": get_surface_type_name(surface_type)
        }
        
        # Extract type-specific parameters
        if surface_type == GeomAbs_Cylinder:
            cylinder = adaptor.Cylinder()
            face_info["radius"] = cylinder.Radius()
            axis = cylinder.Axis()
            face_info["axis"] = {
                "direction": {
                    "x": axis.Direction().X(),
                    "y": axis.Direction().Y(),
                    "z": axis.Direction().Z()
                },
                "location": {
                    "x": axis.Location().X(),
                    "y": axis.Location().Y(),
                    "z": axis.Location().Z()
                }
            }
        
        elif surface_type == GeomAbs_Cone:
            cone = adaptor.Cone()
            face_info["semi_angle"] = math.degrees(cone.SemiAngle())
            face_info["ref_radius"] = cone.RefRadius()
        
        elif surface_type == GeomAbs_Sphere:
            sphere = adaptor.Sphere()
            face_info["radius"] = sphere.Radius()
            center = sphere.Location()
            face_info["center"] = {
                "x": center.X(),
                "y": center.Y(),
                "z": center.Z()
            }
        
        elif surface_type == GeomAbs_Torus:
            torus = adaptor.Torus()
            face_info["major_radius"] = torus.MajorRadius()
            face_info["minor_radius"] = torus.MinorRadius()
        
        elif surface_type == GeomAbs_Plane:
            plane = adaptor.Plane()
            normal = plane.Axis().Direction()
            face_info["normal"] = {
                "x": normal.X(),
                "y": normal.Y(),
                "z": normal.Z()
            }
        
        # Get face area
        props = GProp_GProps()
        brepgprop.SurfaceProperties(face, props)
        face_info["area"] = props.Mass()
        
        faces.append(face_info)
        face_id += 1
        explorer.Next()
    
    return faces


def analyze_edges(shape):
    """Analyze all edges and extract their geometric properties."""
    edges = []
    explorer = TopExp_Explorer(shape, TopAbs_EDGE)
    edge_id = 0
    
    while explorer.More():
        edge = explorer.Current()
        adaptor = BRepAdaptor_Curve(edge)
        curve_type = adaptor.GetType()
        
        edge_info = {
            "id": edge_id,
            "type": get_curve_type_name(curve_type)
        }
        
        # Extract type-specific parameters
        if curve_type == GeomAbs_Circle:
            circle = adaptor.Circle()
            edge_info["radius"] = circle.Radius()
            center = circle.Location()
            edge_info["center"] = {
                "x": center.X(),
                "y": center.Y(),
                "z": center.Z()
            }
        
        elif curve_type == GeomAbs_Ellipse:
            ellipse = adaptor.Ellipse()
            edge_info["major_radius"] = ellipse.MajorRadius()
            edge_info["minor_radius"] = ellipse.MinorRadius()
        
        # Get edge length
        props = GProp_GProps()
        brepgprop.LinearProperties(edge, props)
        edge_info["length"] = props.Mass()
        
        edges.append(edge_info)
        edge_id += 1
        explorer.Next()
    
    return edges


def get_face_type_summary(faces):
    """Summarize face types and their counts."""
    summary = {}
    for face in faces:
        face_type = face["type"]
        if face_type not in summary:
            summary[face_type] = {"count": 0, "total_area": 0}
        summary[face_type]["count"] += 1
        summary[face_type]["total_area"] += face["area"]
    return summary


def detect_primitive_shape(faces, topology, mass_props, bbox):
    """Try to detect if the shape is a common primitive."""
    face_summary = get_face_type_summary(faces)
    
    detected = {"type": "complex", "confidence": "low"}
    
    # Detect cylinder: 1 cylindrical face + 2 circular planar faces
    if (face_summary.get("cylinder", {}).get("count") == 1 and 
        face_summary.get("plane", {}).get("count") == 2):
        
        cyl_face = next(f for f in faces if f["type"] == "cylinder")
        radius = cyl_face["radius"]
        
        # Calculate height from bounding box along cylinder axis
        axis = cyl_face.get("axis", {}).get("direction", {})
        if abs(axis.get("z", 0)) > 0.9:  # Vertical cylinder
            height = bbox["dimensions"]["height"]
        elif abs(axis.get("x", 0)) > 0.9:
            height = bbox["dimensions"]["width"]
        else:
            height = bbox["dimensions"]["depth"]
        
        detected = {
            "type": "cylinder",
            "confidence": "high",
            "parameters": {
                "radius": radius,
                "diameter": radius * 2,
                "height": height
            }
        }
    
    # Detect sphere: single spherical face
    elif (face_summary.get("sphere", {}).get("count") == 1 and 
          topology["faces"] == 1):
        
        sphere_face = next(f for f in faces if f["type"] == "sphere")
        detected = {
            "type": "sphere",
            "confidence": "high",
            "parameters": {
                "radius": sphere_face["radius"],
                "diameter": sphere_face["radius"] * 2,
                "center": sphere_face["center"]
            }
        }
    
    # Detect box: 6 planar faces
    elif (face_summary.get("plane", {}).get("count") == 6 and 
          topology["faces"] == 6):
        
        dims = bbox["dimensions"]
        detected = {
            "type": "box",
            "confidence": "high",
            "parameters": {
                "width": dims["width"],
                "depth": dims["depth"],
                "height": dims["height"]
            }
        }
    
    # Detect cone
    elif (face_summary.get("cone", {}).get("count") == 1 and
          face_summary.get("plane", {}).get("count") in [1, 2]):
        
        cone_face = next(f for f in faces if f["type"] == "cone")
        detected = {
            "type": "cone",
            "confidence": "medium",
            "parameters": {
                "semi_angle_degrees": cone_face["semi_angle"],
                "ref_radius": cone_face["ref_radius"]
            }
        }
    
    return detected


def analyze_step_file(filepath, include_details=True):
    """
    Main function to analyze a STEP file and extract all available information.
    
    Args:
        filepath: Path to the STEP file
        include_details: If True, include detailed face/edge analysis
    
    Returns:
        Dictionary containing all extracted information
    """
    print(f"Loading STEP file: {filepath}")
    shape = load_step(filepath)
    
    print("Extracting bounding box...")
    bbox = get_bounding_box(shape)
    
    print("Calculating mass properties...")
    mass_props = get_mass_properties(shape)
    
    print("Counting topology...")
    topology = get_topology_counts(shape)
    
    print("Analyzing faces...")
    faces = analyze_faces(shape)
    face_summary = get_face_type_summary(faces)
    
    print("Detecting primitive shape...")
    primitive = detect_primitive_shape(faces, topology, mass_props, bbox)
    
    result = {
        "file": filepath,
        "bounding_box": bbox,
        "mass_properties": mass_props,
        "topology": topology,
        "face_type_summary": face_summary,
        "detected_primitive": primitive
    }
    
    if include_details:
        print("Analyzing edges...")
        edges = analyze_edges(shape)
        result["faces"] = faces
        result["edges"] = edges
    
    return result


def print_summary(analysis):
    """Print a human-readable summary of the analysis."""
    print("\n" + "=" * 60)
    print("STEP FILE ANALYSIS SUMMARY")
    print("=" * 60)
    
    print(f"\nFile: {analysis['file']}")
    
    # Detected shape
    prim = analysis["detected_primitive"]
    print(f"\nDetected Shape: {prim['type'].upper()} (confidence: {prim['confidence']})")
    if "parameters" in prim:
        print("Parameters:")
        for key, value in prim["parameters"].items():
            if isinstance(value, dict):
                print(f"  {key}: {value}")
            else:
                print(f"  {key}: {value:.4f}")
    
    # Bounding box
    dims = analysis["bounding_box"]["dimensions"]
    print(f"\nBounding Box Dimensions:")
    print(f"  Width:  {dims['width']:.4f}")
    print(f"  Depth:  {dims['depth']:.4f}")
    print(f"  Height: {dims['height']:.4f}")
    
    # Mass properties
    mp = analysis["mass_properties"]
    print(f"\nMass Properties:")
    print(f"  Volume:       {mp['volume']:.4f}")
    print(f"  Surface Area: {mp['surface_area']:.4f}")
    
    # Topology
    topo = analysis["topology"]
    print(f"\nTopology:")
    print(f"  Solids:   {topo['solids']}")
    print(f"  Faces:    {topo['faces']}")
    print(f"  Edges:    {topo['edges']}")
    print(f"  Vertices: {topo['vertices']}")
    
    # Face types
    print(f"\nFace Types:")
    for face_type, info in analysis["face_type_summary"].items():
        print(f"  {face_type}: {info['count']} (total area: {info['total_area']:.4f})")
    
    print("\n" + "=" * 60)


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python step_analyzer.py <path_to_step_file> [--json]")
        print("\nOptions:")
        print("  --json    Output results as JSON instead of summary")
        sys.exit(1)
    
    filepath = sys.argv[1]
    output_json = "--json" in sys.argv
    
    try:
        analysis = analyze_step_file(filepath)
        
        if output_json:
            print(json.dumps(analysis, indent=2))
        else:
            print_summary(analysis)
            
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
