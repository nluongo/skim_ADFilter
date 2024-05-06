# skim_delphes

Macro skim_Delphes.C reads any Delphes (.root) file, looks into Delphes Ttree, and saves TLorentzVectors for jets, b-tagged jets, electrons, muons, photons into a smaller file skimmed_delphes.root.

To run this, use the following command while replacing "input.root", "skimmed_delphes" (output file) by your own filenames:

root -l -b 'skim_Delphes.C("input_delphes.root", "skimmed_delphes.root", Center_of_mass_energy_value)' 

If Center_of_mass_energy_value=13000 GeV, use : 

root -l -b 'skim_Delphes.C("input_delphes.root", "skimmed_delphes.root", 13000)' 

# skim PHYSLITE (ATLAS) files

Macro skim_Delphes_PHYSLITE.C reads any PHYSLITE (.root) file of ATLAS, looks into "Nominal" Ttree, and saves TLorentzVectors for jets, b-tagged jets, electrons, muons, photons into a smaller file skimmed_delphes.root.

Use this command :

root -l -b 'skim_Delphes_PHYSLITE.C("input_PHYSLITE.root", "skimmed_PHYSLITE.root", Center_of_mass_energy_value)' 

Example: 

root -l -b 'skim_Delphes_PHYSLITE.C("PHYSLITE1.root", "out_phys_light.root", 13000)'


# Read skimmed_delphes file 

Macro read_skimmed_delphes.C reads the file skimmed_delphes.root and prints all the TLorentzVectors for all events. 

To run this, use the following command while replacing "skimmed_delphes.root" by your own file name :

root -l -b 'read_skimmed_delphes.C("skimmed_delphes.root")'
