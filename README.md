SAB re-architected ES's code to streamline how it is used to:
1)  use dm1 & dm2 to correct for small WFE's injected by drift in the system (sab_wfc_loops.m)
2)  evaluate contrast for pre- and post-wfc (sab_eval_dm_psf_loops.m)
Both of these functions call pr2_hex_vvc.m to set up the model and the wavelegths they use.

pr2_hex_vvc.m used to be the main script that could do everything, including digging a dark hole.
It is 1000 lines long and has a lot of "if 0" and "return" spread throughtout it,
so it is not meant to be run end-to-end.  However, it also sets
up the coronagraph model for both in- and out-of-band computations and is still used for this purpose (by SAB).

pr2_hex_vvc.m calls a small script written by SAB, setup_paths.m, which tell the code where falco is
and where the data is that is used to set up the models.  It then calls sab2_setup_hex.m, which is 300 lines,
that explicitly sets up the model (i.e. mp structure).  Depending on the intended use (wfc or eval), the wavelength 
band is set to 1500 nm monochromatic or 500 nm with 12 bands.  It then calls falco_flesh_out_workspace(mp).

ES used to call pr2_hex_vvc.m first, and then either sub_wfc_loops.m to do the WFC and sub_eval_dm_psf_loops.m 
to compute the contrast.  SAB turned it around since he is only performing these latter two studies.
