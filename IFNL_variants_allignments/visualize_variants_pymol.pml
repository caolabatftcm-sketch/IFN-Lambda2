# =============================================================================
# PyMOL Script: IFN-lambda Variant Visualization (Human vs Mouse)
# =============================================================================
# Variants mapped:
#   IFNL2: Human R154H  -->  Mouse R147
#   IFNL3: Human H132R  -->  Mouse H124
#   IFNL3: Human E175K  -->  Mouse E168
#
# UniProt / AlphaFold IDs used:
#   Human IFNL2: Q8IZJ0   Mouse IFNL2: Q4VK74
#   Human IFNL3: Q8IZI9   Mouse IFNL3: Q8CGK6
#
# Requirements: PyMOL 2.4+ (for AlphaFold fetch support)
# Run: pymol visualize_variants_pymol.pml
#      OR open PyMOL, then: File > Run Script
# =============================================================================

# --- 0. Clean up any previous session ---
reinitialize

# =============================================================================
# 1. FETCH ALPHAFOLD STRUCTURES
# =============================================================================
# If fetch fails (no internet), comment out and use local PDB files instead:
#   load /path/to/AF-Q8IZJ0-F1-model_v4.pdb, human_IFNL2
#   load /path/to/AF-Q4VK74-F1-model_v4.pdb, mouse_IFNL2
#   load /path/to/AF-Q8IZI9-F1-model_v4.pdb, human_IFNL3
#   load /path/to/AF-Q8CGK6-F1-model_v4.pdb, mouse_IFNL3

fetch Q8IZJ0, human_IFNL2, type=alphafold
fetch Q4VK74, mouse_IFNL2, type=alphafold
fetch Q8IZI9, human_IFNL3, type=alphafold
fetch Q8CGK6, mouse_IFNL3, type=alphafold

# =============================================================================
# 2. INITIAL DISPLAY SETUP
# =============================================================================
hide everything, all
bg_color white
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set sphere_scale, 0.6
set label_size, 14
set label_font_id, 7
set label_color, black

# Show cartoons for all structures
show cartoon, all

# =============================================================================
# 3. COLOR BY pLDDT (AlphaFold confidence score stored in B-factor column)
#    pLDDT > 90  : very high confidence  -> deep blue
#    pLDDT 70-90 : high confidence       -> cyan/teal
#    pLDDT 50-70 : low confidence        -> yellow
#    pLDDT < 50  : very low confidence   -> orange/red
# =============================================================================
spectrum b, blue_cyan_yellow_red, human_IFNL2, minimum=50, maximum=100
spectrum b, blue_cyan_yellow_red, mouse_IFNL2,  minimum=50, maximum=100
spectrum b, blue_cyan_yellow_red, human_IFNL3, minimum=50, maximum=100
spectrum b, blue_cyan_yellow_red, mouse_IFNL3,  minimum=50, maximum=100

# =============================================================================
# 4. PAIRWISE STRUCTURAL ALIGNMENT (mouse onto human)
# =============================================================================
align mouse_IFNL2, human_IFNL2
align mouse_IFNL3, human_IFNL3

# =============================================================================
# 5. DEFINE VARIANT RESIDUE SELECTIONS
# =============================================================================
# IFNL2 variant: R154H (human) / R147 (mouse)
select var_h_IFNL2_R154, human_IFNL2 and resi 154
select var_m_IFNL2_R147, mouse_IFNL2 and resi 147

# IFNL3 variant: H132R (human) / H124 (mouse)
select var_h_IFNL3_H132, human_IFNL3 and resi 132
select var_m_IFNL3_H124, mouse_IFNL3 and resi 124

# IFNL3 variant: E175K (human) / E168 (mouse)
select var_h_IFNL3_E175, human_IFNL3 and resi 175
select var_m_IFNL3_E168, mouse_IFNL3 and resi 168

# Convenience group: all human variant sites
select all_human_variants, var_h_IFNL2_R154 or var_h_IFNL3_H132 or var_h_IFNL3_E175

# Convenience group: all mouse equivalent sites
select all_mouse_variants, var_m_IFNL2_R147 or var_m_IFNL3_H124 or var_m_IFNL3_E168

# =============================================================================
# 6. VISUALIZE VARIANT RESIDUES
# =============================================================================
# Show as sticks (so side chain is visible) + sphere on Calpha
show sticks,  all_human_variants or all_mouse_variants
show spheres, all_human_variants or all_mouse_variants
set stick_radius, 0.15
set sphere_scale, 0.35, all_human_variants or all_mouse_variants

# Color coding:
#   Human variants: warm colors   Mouse equivalents: cool colors
color firebrick,  var_h_IFNL2_R154    # Human IFNL2 R154 - dark red
color salmon,     var_m_IFNL2_R147    # Mouse IFNL2 R147 - light red/salmon

color marine,     var_h_IFNL3_H132    # Human IFNL3 H132 - dark blue
color lightblue,  var_m_IFNL3_H124    # Mouse IFNL3 H124 - light blue

color forest,     var_h_IFNL3_E175    # Human IFNL3 E175 - dark green
color palegreen,  var_m_IFNL3_E168    # Mouse IFNL3 E168 - light green

# =============================================================================
# 7. LABELS
# =============================================================================
label var_h_IFNL2_R154 and name CA, "R154 (Hs IFNL2)"
label var_m_IFNL2_R147 and name CA, "R147 (Mm IFNL2)"
label var_h_IFNL3_H132 and name CA, "H132 (Hs IFNL3)"
label var_m_IFNL3_H124 and name CA, "H124 (Mm IFNL3)"
label var_h_IFNL3_E175 and name CA, "E175 (Hs IFNL3)"
label var_m_IFNL3_E168 and name CA, "E168 (Mm IFNL3)"

# =============================================================================
# 8. PRINT pLDDT VALUES AT VARIANT POSITIONS (runs in PyMOL Python console)
# =============================================================================
python
cmd.iterate("var_h_IFNL2_R154 and name CA", "print(f'Hs IFNL2 R154  pLDDT={b:.1f}')")
cmd.iterate("var_m_IFNL2_R147 and name CA", "print(f'Mm IFNL2 R147  pLDDT={b:.1f}')")
cmd.iterate("var_h_IFNL3_H132 and name CA", "print(f'Hs IFNL3 H132  pLDDT={b:.1f}')")
cmd.iterate("var_m_IFNL3_H124 and name CA", "print(f'Mm IFNL3 H124  pLDDT={b:.1f}')")
cmd.iterate("var_h_IFNL3_E175 and name CA", "print(f'Hs IFNL3 E175  pLDDT={b:.1f}')")
cmd.iterate("var_m_IFNL3_E168 and name CA", "print(f'Mm IFNL3 E168  pLDDT={b:.1f}')")
python end

# =============================================================================
# 9. SAVE SCENES FOR EACH PROTEIN PAIR
# =============================================================================
# Scene 1: IFNL2 pair overview
disable human_IFNL3
disable mouse_IFNL3
orient human_IFNL2 or mouse_IFNL2
scene IFNL2_pair, store

# Scene 2: IFNL2 variant close-up
zoom var_h_IFNL2_R154 or var_m_IFNL2_R147, 8
scene IFNL2_R154_closeup, store

# Scene 3: IFNL3 pair overview
enable human_IFNL3
enable mouse_IFNL3
disable human_IFNL2
disable mouse_IFNL2
orient human_IFNL3 or mouse_IFNL3
scene IFNL3_pair, store

# Scene 4: IFNL3 H132/H124 close-up
zoom var_h_IFNL3_H132 or var_m_IFNL3_H124, 8
scene IFNL3_H132_closeup, store

# Scene 5: IFNL3 E175/E168 close-up
zoom var_h_IFNL3_E175 or var_m_IFNL3_E168, 8
scene IFNL3_E175_closeup, store

# Restore all
enable all
orient all
scene all_structures, store

# =============================================================================
# 10. EXPORT IMAGES (optional — comment out if not needed)
# =============================================================================
# Uncomment below to auto-export PNG images at 300 DPI (2400x2400 px)
# set ray_opaque_background, on
# ray 2400, 2400
# png IFNL2_pair.png, dpi=300

# =============================================================================
# DONE — Use Scene panel (top menu: Display > Scenes) to toggle views.
# Useful PyMOL commands to explore further:
#   set_view / get_view        - save/restore camera positions
#   dssp                       - assign secondary structure
#   show surface, human_IFNL2  - surface representation
#   color white, human_IFNL2   - recolor before applying surface
# =============================================================================
print "=== IFN-lambda variant script loaded successfully ==="
print "Scenes saved: IFNL2_pair, IFNL2_R154_closeup, IFNL3_pair, IFNL3_H132_closeup, IFNL3_E175_closeup, all_structures"
print "pLDDT values printed above for each variant position."
