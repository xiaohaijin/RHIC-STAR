(TeX-add-style-hook
 "STARBranch"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("beamer" "10pt")))
   (add-to-list 'LaTeX-verbatim-environments-local "semiverbatim")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "beamer"
    "beamer10"
    "ctex")
   (TeX-add-symbols
    "meeting")
   (LaTeX-add-labels
    "BTofHeader"
    "BTofHit"
    "BTofRawHit"
    "CovGlobTrack"
    "CovPrimTrack"
    "CpvCluster"
    "CpvHit"
    "DetectorStates"
    "EEmcPrs"
    "EEmcSmdu"
    "EEmcSmdv"
    "EmcPrs"
    "EmcSmde"
    "EmcSmdp"
    "EmcTow"
    "Event"
    "FmsHit"
    "GlobalTracks"
    "KinkAssoc"
    "Kink"
    "L3AlgoAccept"
    "L3AlgoReject"
    "L3Tracks"
    "McEvent"
    "McKink"
    "McV0"
    "McXi"
    "MuEvent"
    "OtherTracks"
    "PmdCluster"
    "PmdHit"
    "pp2pp"
    "PrimaryTracks"
    "PrimaryVertices"
    "RichSpectra"
    "StrangeCuts"
    "TofData"
    "TofRawData"
    "V0Assoc"
    "V0"
    "XiAssoc"
    "Xi")
   (LaTeX-add-xcolor-definecolors
    "seagreen"))
 :latex)

