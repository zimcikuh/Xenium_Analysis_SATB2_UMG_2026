samples_map = {
    "C1": "Satb2_+/+_1",
    "C2": "Satb2_+/+_2",
    "C3": "Satb2_+/+_3",
    "C4": "Satb2_+/+_4",
    "C5": "Satb2_-/-_1",
    "C6": "Satb2_-/-_2",
    "C7": "Satb2_-/-_3",
    "C8": "Satb2_-/-_4",
}

data_markers = [
    {
        "final_annotation": "Craniofacial developmental mesenchyme",
        "final_annotation_code": "craniofacial_developmental_mesenchyme",
        "lineage": "Mesenchymal",
        "key_markers": ["Pax3", "Osr1", "Lhx8", "Cxcl12", "Wnt5a"],
        "description": "Embryonic craniofacial tooth patterning rather",
    },
    {
        "final_annotation": "Osteogenic progenitors",
        "final_annotation_code": "osteogenic_progenitors",
        "lineage": "Osteogenic",
        "key_markers": ["Msx1", "Runx2", "Alx1", "Fzd1", "Thbs2"],
        "description": "Mesenchymal cells with a progenitor state upstr",
    },
    {
        "final_annotation": "Myogenic progenitors",
        "final_annotation_code": "myogenic_progenitors",
        "lineage": "Muscle",
        "key_markers": ["Pax7", "Fgfr4", "Myod1", "Myog", "Chrna1"],
        "description": "Skeletal muscle-lineage captured in the ROI.",
    },
    {
        "final_annotation": "Differentiated skeletal muscle cells",
        "final_annotation_code": "differentiated_skeletal_muscle_cells",
        "lineage": "Muscle",
        "key_markers": ["Cacna1s", "Trdn", "Casq2", "Mybph", "Chrng"],
        "description": "Contractile skeletal mus surrounding tissue.",
    },
    {
        "final_annotation": "ECM-producing mesenchymal cells",
        "final_annotation_code": "ECM_producing_mesenchymal_cells",
        "lineage": "Mesenchymal",
        "key_markers": ["Tnxb", "Mfap5", "Barx1", "Sfrp2", "Igf1", "Col14a1"],
        "description": "Mesenchymal cells enric structural support and si",
    },
    {
        "final_annotation": "Mesenchymal progenitors",
        "final_annotation_code": "mesenchymal_progenitors",
        "lineage": "Mesenchymal",
        "key_markers": ["Twist1", "Twist2", "Tbx18", "Ntn1", "Tgfbr3"],
        "description": "Undifferentiated mesenc stromal pool that can fee",
    },
    {
        "final_annotation": "Basal epithelial cells",
        "final_annotation_code": "basal_epithelial_cells",
        "lineage": "Epithelial",
        "key_markers": ["Lamb3", "Lama3", "Dsg2", "Itgb4", "Dsc3", "Grhl2"],
        "description": "Basal epithelial cells with structural epithelial comp",
    },
    {
        "final_annotation": "Endothelial cells",
        "final_annotation_code": "endothelial_cells",
        "lineage": "Vascular",
        "key_markers": ["Aplnr", "Nos3", "Esam", "Flt1", "Pecam1", "Cdh5"],
        "description": "Canonical vascular endot main blood-vessel compa",
    },
    {
        "final_annotation": "Differentiated osteoblast-lineage cells",
        "final_annotation_code": "differentiated_osteoblast_lineage_cells",
        "lineage": "Osteogenic",
        "key_markers": ["Sp7", "Ibsp", "Pth1r", "Smpd3", "Col24a1"],
        "description": "Differentiating osteoblast developing bone-forming",
    },
    {
        "final_annotation": "Keratinized epithelial cells",
        "final_annotation_code": "keratinized_epithelial_cells",
        "lineage": "Epithelial",
        "key_markers": ["Acer1", "Dsc1", "Pla2g2f", "Pnpla1", "Lipm"],
        "description": "Highly differentiated kerat mature epithelial surface",
    },
    {
        "final_annotation": "Early odontogenic mesenchyme",
        "final_annotation_code": "early_odontogenic_mesenchyme",
        "lineage": "Odontogenic",
        "key_markers": ["Msx1", "Msx2", "Pax9", "Dlx2", "Dlx3", "Bmp4", "Fgf10", "Runx2"],
        "description": "Tooth-competent mesenc without mature osteoblast patterning.",
    },
    {
        "final_annotation": "Dental epithelial cells",
        "final_annotation_code": "dental_epithelial_cells",
        "lineage": "Epithelial",
        "key_markers": ["Pitx2", "Dlx1", "Dlx2", "Msx2", "Wnt7a"],
        "description": "Dental epithelium with pat participates in epithelial cc",
    },
    {
        "final_annotation": "Matrix-remodeling osteogenic cells",
        "final_annotation_code": "matrix_remodeling_osteogenic_cells",
        "lineage": "Osteogenic",
        "key_markers": ["Mmp13", "Mmp9", "Dmp1", "Enpp1", "Fam20c", "Ank"],
        "description": "Osteogenic cells with stron active remodeling state rat",
    },
    {
        "final_annotation": "SATB2-positive signaling mesenchymal cells",
        "final_annotation_code": "SATB2_positive_signaling_mesenchymal_cells",
        "lineage": "Mesenchymal",
        "key_markers": ["Sfrp4", "Rspo1", "Inhba", "Alx1", "Satb2"],
        "description": "Signaling-active mesenchym regulatory mesenchymal st",
    },
    {
        "final_annotation": "Signaling epithelial cells",
        "final_annotation_code": "signaling_epithelial_cells",
        "lineage": "Epithelial",
        "key_markers": ["Tfap2b", "Cxcl14", "Wnt6", "Fzd10", "Gata3", "Pou2f3"],
        "description": "Epithelial cells enriched for regulatory epithelial state v",
    },
    {
        "final_annotation": "Tendon/ligament-like mesenchymal cells",
        "final_annotation_code": "tendon_ligament_like_mesenchymal_cells",
        "lineage": "Mesenchymal",
        "key_markers": ["Scx", "Mkx", "Thbs4", "Fmod", "Meox1", "Meox2"],
        "description": "Mesenchymal cells with ten ligament-like or other fibrou",
    },
    {
        "final_annotation": "Macrophages",
        "final_annotation_code": "macrophages",
        "lineage": "Immune",
        "key_markers": ["Fcrls", "Mrc1", "Adgre1", "Csf1r", "Cx3cr1", "Aif1"],
        "description": "Myeloid/macrophage popul contributes to local phagocy",
    },
    {
        "final_annotation": "Mature osteoblast-lineage cells",
        "final_annotation_code": "mature_osteoblast_lineage_cells",
        "lineage": "Osteogenic",
        "key_markers": ["Bglap", "Bglap2", "Dmp1", "Phex", "Ibsp"],
        "description": "More mature osteoblast-lin Represents a later stage tha",
    },
    {
        "final_annotation": "Differentiated epithelial cells",
        "final_annotation_code": "differentiated_epithelial_cells",
        "lineage": "Epithelial",
        "key_markers": ["Elf3", "Gjb2", "Hopx", "Csta1", "Serpinb3b"],
        "description": "Differentiated epithelial cell still consistent with a matur",
    },
    {
        "final_annotation": "Perivascular smooth muscle cells",
        "final_annotation_code": "perivascular_smooth_muscle_cells",
        "lineage": "Vascular",
        "key_markers": ["Myocd", "Actg2", "Osr2", "Tbx18", "Cdh6"],
        "description": "Smooth-muscle-like mural c and local stromal interaction",
    },
    {
        "final_annotation": "Neural crest-derived cells",
        "final_annotation_code": "neural_crest_derived_cells",
        "lineage": "Neural crest",
        "key_markers": ["Sox10", "Foxd3", "Ednrb", "L1cam", "Gfra3"],
        "description": "Neural crest-derived / Schwa specialized non-mesenchym",
    },
    {
        "final_annotation": "Vascular smooth muscle cells",
        "final_annotation_code": "vascular_smooth_muscle_cells",
        "lineage": "Vascular",
        "key_markers": ["Myh11", "Actg2", "Kcnj8", "Abcc9", "Gucy1a1", "Tagln"],
        "description": "Canonical contractile vascula vascular smooth muscle com",
    },
    {
        "final_annotation": "Lymphatic endothelial cells",
        "final_annotation_code": "lymphatic_endothelial_cells",
        "lineage": "Vascular",
        "key_markers": ["Ccl21a", "Ccl21b", "Lyve1", "Flt4", "Prox1"],
        "description": "Lymphatic endothelial cells v Distinct from the blood endo",
    },
    {
        "final_annotation": "Mast cells",
        "final_annotation_code": "mast_cells",
        "lineage": "Immune",
        "key_markers": ["Tpsb2", "Hdc", "Il1rl1", "Gata2", "Kit"],
        "description": "Mast cell population with his immune compartment.",
    },
    {
        "final_annotation": "Chondrocytes",
        "final_annotation_code": "chondrocytes",
        "lineage": "Chondrogenic",
        "key_markers": ["Col9a2", "Col9a3", "Acan", "Matn1", "Hapln1"],
        "description": "Cartilage-forming cells with a cartilaginous developmental",
    },
    {
        "final_annotation": "Osteoclasts",
        "final_annotation_code": "osteoclasts",
        "lineage": "Immune/Osteoclast ic",
        "key_markers": ["Acp5", "Atp6v0d2", "Ctsk", "Mmp9", "Tnfrsf11a"],
        "description": "Bone-resorbing osteoclasts v macrophages despite shared",
    },
    {
        "final_annotation": "Epithelial cells",
        "final_annotation_code": "epithelial_cells",
        "lineage": "Epithelial",
        "key_markers": ["Epcam", "Cldn7", "Elf3", "Il33"],
        "description": "Small epithelial population w if finer epithelial subclustering",
    },
    {
        "final_annotation": "Angiogenic endothelial cells",
        "final_annotation_code": "angiogenic_endothelial_cells",
        "lineage": "Vascular",
        "key_markers": ["Kdr", "Cdh5", "Dll4", "Notch4", "Cd93", "Podxl"],
        "description": "Activated endothelial cells wi endothelial state rather than",
    },
    {
        "final_annotation": "Hypertrophic chondrocytes",
        "final_annotation_code": "hypertrophic_chondrocytes",
        "lineage": "Chondrogenic",
        "key_markers": ["Col10a1", "Acan", "Mmp13", "Fam20c", "Ccn2"],
        "description": "Hypertrophic chondrocytes at features. Likely bridge cartilag",
    },
]




# Keys match the values of cell_type_mapping EXACTLY (case + spacing matter)
PALETTE = {
  # ── Epithelial (greens) ──────────────────────────────────────────────────
  'Epithelial': {
    'Basal epithelial cells':                          '#0B4F2A',   # deep forest
    'Oral epithelial cells':                           '#2DB84C',   # vivid kelly
    'Glandular epithelial cells':                      '#558B2F',   # moss '#558B2F',
    'Dental epithelial cells':                         "#8FBE62",   # olive   '#7A7A1F', 
    'Signaling epithelial cells':                      '#7A7A1F',   # light grass A6E36B
    'Keratinized epithelial cells':                    '#9DB39A',   # sage
    'Differentiated epithelial cells':                 "#D1E9B6"},  # pale green
  # ── Mesenchymal (purples / magenta) ──────────────────────────────────────
  'Mesenchymal': {
    'Mesenchymal progenitor cells':                    '#5B21B6',   # deep violet
    'Craniofacial developmental mesenchymal cells':    '#8B5CF6',   # violet
    'ECM producing mesenchymal cells':                 '#C026D3',   # magenta
    'Outer fibroperiosteal cells':                     '#C4B5FD'},  # lavender
  # ── Dental mesenchyme (golds; follicle -> papilla) ───────────────────────
  'Dental mesenchyme': {
    'Dental follicle mesenchymal cells':               '#8C6D1F',   # bronze
    'Dental follicle cells':                           '#C9A227',   # old gold
    'Dental follicle signaling niche':                 '#F2C94C',   # gold
    'Dental papilla cells':                            '#FFE699',   # pale butter
    'Mut-specific dental cells':                       '#1A1A1A'},  # near-black – highlight
  # ── Neural crest (warm grey) ─────────────────────────────────────────────
  'Neural crest': {
    'Neural crest':                                    '#78716C'},
  # ── Osteogenic (blues; progenitor light -> mature navy) ──────────────────
  'Osteogenic': {
    'Osteogenic progenitors':                          '#22D3EE',   # cyan
    'Early osteoblasts':                               '#14B8A6',   # teal
    'Matrix-secreting osteoblast lineage cells':       '#0F766E',   # dark teal
    'Pre-osteocytes / diff. osteoblast lineage cells': '#134E4A'},  # deep teal
  # ── Muscle (browns) ──────────────────────────────────────────────────────
  'Muscle': {
    'Myogenic progenitor':                             '#C08A5B',   # tan
    'Skeletal muscle':                                 '#6D3B1E'},  # dark brown
  # ── Chondrogenic (reds) ──────────────────────────────────────────────────
  'Chondrogenic': {
    'Chondrocytes':                                    '#FF1744',   # vivid scarlet
    'Hypertrophic chondrocytes':                       '#B71C1C'},  # deep red
  # ── Vascular / mural (teals; EC light -> SMC dark) ───────────────────────
  'Vascular': {
    'Endothelial cells':                               '#42A5F5',   # mid blue
    'Lymphatic endothelial cells':                     '#8FD3FF',   # light sky
    'perivascular smooth muscle cells':                '#084EFF',   # vivid royal
    'Vascular smooth muscle cells':                    '#0B1F5C'},  # deep navy
  # ── Immune (oranges + mast rose) ─────────────────────────────────────────
  'Immune': {
    'Macrophages':                                     '#FF8A3D',   # orange
    'Osteoclasts':                                     '#FF6B9D',   # burnt orange C2410C
    'Mast cells':                                      '#C2410C'},  # rose FF6B9D
}
 
LEVEL1_COLORS = {
    'Epithelial':        '#2DB84C',
    'Mesenchymal':       '#8B5CF6',
    'Dental mesenchyme': '#C9A227',
    'Neural crest':      '#78716C',
    'Osteogenic':        '#14B8A6',
    'Muscle':            '#6D3B1E',
    'Chondrogenic':      '#FF1744',
    'Vascular':          '#0B1F5C',
    'Immune':            '#FF8A3D',
}
 
# flat lookups
L2_COLORS = {t: c for g in PALETTE.values() for t, c in g.items()}
L2_TO_L1  = {t: g for g, d in PALETTE.items() for t in d}
L2_ORDER  = list(L2_COLORS)   # palette order = grouped by lineage
 