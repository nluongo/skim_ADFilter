# skim_delphes

Macro skim_Delphes.C reads any Delphes (.root) file, looks into Delphes Ttree, and saves TLorentzVectors for jets, b-tagged jets, electrons, muons, photons into a smaller file skimmed_delphes.root.

To run this, use the following command while replacing "input.root", "skimmed_delphes" (output file) by your own filenames:

root -l -b 'skim_Delphes.C("input.root", "skimmed_delphes.root")'

# Read skimmed_delphes file 

Macro read_skimmed_delphes.C reads the file skimmed_delphes.root and prints all the TLorentzVectors for all events. 

To run this, use the following command while replacing "skimmed_delphes.root" by your own file name :

root -l -b 'read_skimmed_delphes.C("skimmed_delphes.root")'
