# heavify

Calculate a heavy neutrino flux from a dk2nu format neutrino beamline simulation.

This creates a TH2D named `h_HN{mu/e}_d{DD}_{LLL}`, with neutrino momentum on the x-axis, and heavy neutrino mass on the y-axis. `dDD` is the decay type (ndecay in dk2nu format, e.g. `d05` is K+ decay). `LLL` is a location, defined in the dkMeta tree of the input.

The normalisation is such that the coupling constant Θ is 1. A different flux is created for the muon coupling `h_HNmu_...` and electron coupling `h_HNe_...`, to allow your analysis to change the coupling constants independently.

Also creates a TH1D of the light neutrino flux, for comparison.

Also creates a TH2D named `h_ndecay_nu{mu/e}` of the decay types contributing to the heavy neutrino creation

## Compiling

Either compile directly, linking to ROOT and DK2NU.

Or compile within root: `root -l -b heavify.C`
Trying to compile directly, `.L heavify.C+`, will cause errors, you need to define the DK2NU include paths first.
The macro assume you have loaded DK2NU via ups first. If using a standalone DK2NU, modify the last section of the program from the line `#else // if __CINT__ is defined`.
