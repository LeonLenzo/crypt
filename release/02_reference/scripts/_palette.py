"""The one palette for the two panels. figures/tree/tree.R and figures/tree/tree_legend.R
mirror the MATRIX hexes in their CAT_COLOUR -- keep all three in step.

ONE scheme across both panels, ordered broad -> specific, cold -> warm. The scope funnel runs
dark blue (every plant record) to light blue (has a genome); the buildable set is grey (the
hinge: it is the database, and the database is mostly residents); then the label categories
warm up green -> yellow -> orange as they get more specific, ending at the wheat pathogen.
Panel A and panel B share grey (buildable = resident) and orange (its speck = the tree's
wheat-pathogen tips), so the two panels read as one figure.
"""

# the collapse stages, broadest (darkest) to buildable
PLANT_RECORD = "#08306b"   # recorded on a plant (USDA)      -- darkest blue
KNOWN_TAXON  = "#2171b5"   # a known NCBI taxon              -- dark blue
HAS_GENOME   = "#9ecae1"   # has a genome assembly           -- light blue

# the four label categories (and the two panel A shares)
RESIDENT   = "#d9d9d9"     # annotated buildable / resident  -- grey   (both panels)
ON_WHEAT   = "#41ab5d"     # recorded on wheat, not pathogen -- green
PLANT_PATH = "#fdd130"     # known plant pathogen, not wheat -- yellow
WHEAT_PATH = "#e8702a"     # known wheat pathogen (PHI-base) -- orange (both panels)

INK   = "#1a1a1a"
MUTED = "#707070"


# ── the matrix (panel A) ──────────────────────────────────────────────────────
# A 2-D ramp, because the matrix crosses two independent axes and colour has to carry both.
# HUE is the PHI-base column (green non-pathogen -> amber plant pathogen -> orange wheat
# pathogen); SATURATION is the USDA row (saturated = recorded on wheat, pale = recorded on
# some other plant, grey = congener). That retires the fill-plus-outline encoding and the
# legend with it: a cell's position already states both facts, and its colour agrees.
MATRIX = {
    ("wheat", "none"):       "#1aa64b",
    ("wheat", "vagrant"):    "#f5a800",
    ("wheat", "wheat_path"): "#e05a26",
    ("plant", "none"):       "#c3e8bf",
    ("plant", "vagrant"):    "#fbe49a",
    ("congener", "none"):    "#f0f0f0",
}


def darken(hex_colour, f=0.55):
    """Cell borders are the fill walked toward black, so each keeps its own hue family."""
    h = hex_colour.lstrip("#")
    return "#%02x%02x%02x" % tuple(round(int(h[i:i + 2], 16) * f) for i in (0, 2, 4))
