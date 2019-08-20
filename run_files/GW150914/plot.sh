pycbc_inference_plot_posterior --verbose \
    --input-file ./gw150914_posteriors_thinned.hdf \
    --output-file q_chieff_gw150914.png \
    --plot-density \
    --density-cmap Purples \
    --z-arg logplr \
    --plot-contours \
    --contour-color purple \
    --plot-marginal \
    --parameters "(primary_mass(mass1, mass2))/(secondary_mass(mass1, mass2)):\$q$" \
                 chi_eff \
    --contour-color "purple" \
    --maxs "(primary_mass(mass1, mass2))/(secondary_mass(mass1, mass2)):2"



pycbc_inference_plot_posterior --verbose \
    --input-file ./gw150914_posteriors_thinned.hdf \
    --output-file iota_dl_gw150914.png \
    --plot-density \
    --density-cmap Purples \
    --z-arg logplr \
    --plot-contours \
    --contour-color purple \
    --plot-marginal \
    --parameters "inclination*180/pi:$\iota$ (deg)" \
                  distance \
    --mins "inclination*180/pi:0" \
    # --maxs "inclination*180/pi:180" \