(TeX-add-style-hook
 "mesh_one_bud_scar"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("standalone" "tikz" "crop" "convert={density=200,outext=.png}" "border=0.0cm")))
   (TeX-run-style-hooks
    "latex2e"
    "standalone"
    "standalone10"
    "amsmath"
    "physics"
    "graphicx"
    "xcolor")
   (LaTeX-add-xcolor-definecolors
    "boundary_hole"
    "hole"
    "adjacent"
    "rest"))
 :latex)

