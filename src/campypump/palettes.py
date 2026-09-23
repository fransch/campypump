"""One palette for every campypump figure.

Base colours:

    #DC143C crimson   #D8BFD8 thistle   #4169E1 royal blue   #72B272 green   #FFD700 gold

Two colours are derived from the base set because the originals do not work as line or marker
colours: THISTLE_DEEP darkens #D8BFD8, which is too pale to read as a 1.5 px branch, and the two
light variants pair with their parent colour for related categories (chicken / chicken meat,
cattle / beef). No hue carries two meanings within a figure, and only one blue is used.
"""

# --- base ------------------------------------------------------------------
CRIMSON = "#DC143C"
THISTLE = "#D8BFD8"
ROYAL = "#4169E1"
GREEN = "#72B272"
GOLD = "#FFD700"

# --- derived ---------------------------------------------------------------
THISTLE_DEEP = "#8A6FA8"   # darkened thistle, for line work
GOLD_LIGHT = "#FFE97A"
GREEN_LIGHT = "#A8D8A8"
CRIMSON_LIGHT = "#E88BA0"
GREY = "#9A9A9A"
GREY_PALE = "#D9D9D9"
INK = "#111111"

# --- semantic assignments --------------------------------------------------
# cmeB type: the classification itself
TYPE = {
    "RE-cmeB": CRIMSON,
    "wildtype-cmeB": ROYAL,
    "mosaic": GOLD,
    "ambiguous": GREY,
    "unresolved": GREY,
    "not_found": GREY_PALE,
    "incomplete": GREY_PALE,
}
TYPE_SHORT = {"RE": CRIMSON, "WT": ROYAL, "MOS": GOLD, "AMB": GREY, "UNR": GREY,
              "REF": INK, "OUT": THISTLE_DEEP, "mixed": GREY_PALE}

# species: green and deep thistle, so they never collide with the type colours above
SPECIES = {"Cc": GREEN, "Cj": THISTLE_DEEP, "REF": INK, "OUT": GREY, "mixed": GREY_PALE}
SPECIES_LINE = [GREEN, THISTLE_DEEP, CRIMSON, ROYAL]
SPECIES_FILL = {"Cc": GREEN_LIGHT, "Cj": THISTLE}

# isolation source: related categories share a hue at two lightnesses
# THISTLE_DEEP is reserved for C. jejuni, so pigs take the crimson pair instead: within one figure
# no colour may mean two things.
SOURCES = [
    ("chicken offal or meat", GOLD_LIGHT),
    ("chicken", GOLD),
    ("human stool", ROYAL),
    ("pork offal or meat", CRIMSON_LIGHT),
    ("pig", CRIMSON),
    ("beef offal or meat", GREEN_LIGHT),
    ("cattle faeces", GREEN),
    ("wild bird", THISTLE),
]
SOURCE_OTHER = GREY_PALE

# segment plots
SEGMENT_RE = CRIMSON
SEGMENT_WT = ROYAL
