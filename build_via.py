from pyaedt import Hfss
import os

# --- 1. INITIALIZE HFSS ---
# This will open AEDT 2023.x (or your student version)
hfss = Hfss(specified_version="2023.1", non_graphical=False, new_desktop_session=True)

# --- 2. DEFINE VARIABLES ---
# Using variables allows for easy optimization later
hfss["v_drill"] = "0.2mm"
hfss["v_pad"] = "0.4mm"
hfss["v_ap"] = "0.7mm"
hfss["sub_h"] = "1.0mm"
hfss["d_spacing"] = "0.5mm"

# --- 3. CREATE SUBSTRATE ---
# Centering the box around the origin
hfss.modeler.create_box(
    position=["-2mm", "-2mm", "0mm"],
    dimensions_list=["4mm", "4mm", "sub_h"],
    name="Substrate",
    matname="FR4_epoxy"
)

# --- 4. CREATE DIFFERENTIAL VIAS ---
# Positive Via
hfss.modeler.create_cylinder(
    cs_axis="Z",
    position=["d_spacing/2", "0", "0"],
    radius="v_drill/2",
    height="sub_h",
    name="Via_P",
    matname="copper"
)

# Negative Via
hfss.modeler.create_cylinder(
    cs_axis="Z",
    position=["-d_spacing/2", "0", "0"],
    radius="v_drill/2",
    height="sub_h",
    name="Via_N",
    matname="copper"
)

# --- 5. CREATE ANTIPADS (Subtracting from Ground) ---
# Create a Ground Plane
gnd = hfss.modeler.create_rectangle(
    cs_plane="XY",
    position=["-2mm", "-2mm", "0mm"],
    dimension_list=["4mm", "4mm"],
    name="GND_Top"
)

# Create two circles to subtract
ap_p = hfss.modeler.create_circle(cs_plane="XY", position=["d_spacing/2", "0", "0"], radius="v_ap/2")
ap_n = hfss.modeler.create_circle(cs_plane="XY", position=["-d_spacing/2", "0", "0"], radius="v_ap/2")

# Subtract circles from Ground Plane to make the antipad holes
hfss.modeler.subtract(gnd, [ap_p, ap_n], keep_originals=False)

# --- 6. SETUP ANALYSIS ---
setup = hfss.create_setup(setupname="MySetup")
setup.props["Frequency"] = "12GHz"
hfss.create_linear_count_sweep(
    setupname="MySetup",
    unit="GHz",
    freqstart=0.1,
    freqstop=40,
    num_of_points=401,
    sweepname="Sweep1",
    save_fields=False
)

print("Via Model built successfully. Check HFSS window.")