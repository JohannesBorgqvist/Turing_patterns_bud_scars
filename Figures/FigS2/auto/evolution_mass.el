(TeX-add-style-hook
 "evolution_mass"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("standalone" "tikz" "border=0mm")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("fontenc" "T1") ("babel" "swedish" "english")))
   (TeX-run-style-hooks
    "latex2e"
    "standalone"
    "standalone10"
    "fontenc"
    "babel"
    "amsmath"
    "tikz"
    "amssymb"
    "pgfplots")
   (TeX-add-symbols
    "HUGE"
    "hugest"))
 :latex)

