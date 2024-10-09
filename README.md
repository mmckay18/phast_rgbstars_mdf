# M31 Red Giant Branch Stars Metallicity Distribution Function and Gradient

By Myles McKay(@mmckay18),

![Astronomy](https://img.shields.io/badge/Field-Astronomy-blue)
![Last Commit](https://img.shields.io/github/last-commit/mmckay18/phast_rgbstars_mdf)
![Repo Size](https://img.shields.io/github/repo-size/mmckay18/phast_rgbstars_mdf)
![Contributors](https://img.shields.io/github/contributors/mmckay18/phast_rgbstars_mdf)
![Pull Requests](https://img.shields.io/github/issues-pr/mmckay18/phast_rgbstars_mdf)
![License](https://img.shields.io/github/license/mmckay18/phast_rgbstars_mdf)
![Issues](https://img.shields.io/github/issues/mmckay18/phast_rgbstars_mdf)
![Stars](https://img.shields.io/github/stars/mmckay18/phast_rgbstars_mdf)
![Forks](https://img.shields.io/github/forks/mmckay18/phast_rgbstars_mdf)

# Table of Contents

- [Introduction](#introduction)
- [Research Abstract](#research-abstract)
- [Data Sources](#data-sources)
- [Results and Evaluation](#results-and-evaluation)
- [References](#references)
- [License](#license)

# Introduction

- This repo demonstrates the exploratory data analysis, statistical method and viusalization of PHAT and PHAST Red Giant Branch stars.

# Research Abstract

The Panchromatic Hubble Andromeda Southern Treasury (PHAST) (Chen et al. 2024) is an extension of the 2012 Panchromatic Hubble Andromeda Treasury (PHAT) (Dalcanton et al. 2012) that built a rich photometry catalog of the stars in the northern half of M31. We use the optical photomtry of the PHAST to select the RGB stellar population and isochrone models to estimate their stellar metallicity. The stars are then binned to $0.01^\circ$ square bins and we compute the median for metallicity. The median metallicity binned maps sampled in even elliptical bins and compute the median from the inner to outer radii to create the RGB metallicity gradient. Using the stellar photometry of RGB stars in M31, we measure the metallicity gradient in both the northern and southern halves and compare their gradients. The stars comes from the PHAT and PHAST HST surveys optical bands F475W and F814W. Metallicity is calculated from isochrone models using 2D linear interpolation methods, and an RGB selection box is applied to mitigate photometric bias and ensure completeness. Since dust in galaxies affects photometry and can lead to overestimation of metallicity, we use infrared observations to correct for dust effects. While the PHAT catalog includes IR observations, the PHAST catalog only has optical observations, so we use previous IR observations of M31 to measure high dust mass surface density. A previous study of the PHAT metallicity gradient reported a slope of $-0.02 \text{dex/kpc}$ for the northern half of M31, indicating that the RGB population is well-mixed throughout the disk, assuming a fiducial age of $4 Gyr$. In this chapter, we present the metallicity analysis of the PHAST catalog and the comprehensive RGB gradient obtained by combining the PHAT and PHAST catalogs. We find that M31 has a shallow metallicity gradient, with similar slopes $(-0.01 dex/kpc)$ for both the northern and southern halves, suggesting that the RGB stellar population is well-mixed throughout the disk.

# Data

## Data Sources

### Photometry Data

PHAT - https://archive.stsci.edu/hlsp/phat

PHAST - https://www.stsci.edu/cgi-bin/get-proposal-info?observatory=HST&id=16796

### Isochrone Model Table Data

CMD - http://stev.oapd.inaf.it/cgi-bin/cmd

### M31 Dust Map

Drain Spitzer Dust Mass Surface Density Maps - https://www.astro.princeton.edu/~draine/m31dust/m31dust.html

## Data Acquisition

## Data Preprocessing

1. Apply good star criteria as outline in the [2]
2. Combine the photometry catalogs for PHAT and PHAST
3. Apply photometric cut to remove stars dimmer that 23 mag in F814W filter for the isochrone table and the photometry table
4. Use linear interpolation to match isochrone table and photometry table

# Code Structure

# Results and Evaluation

## PHAT and PHAST RGB Spatial Metallicity Map

![Spatial Metallicity Map](/images/m31_analysis_map_subplot.pdf)

## RGB Metallicity Distrubution Function

## RGB Metallicity Gradient

# Future Work

# References

[1] Dalcanton, J. J., Williams, B. F., Lang, D., et al. 2012, The Astrophysical Journal Supplement Series, 200, 18, doi: 10.1088/0067-0049/200/2/18

[2] Williams, B. F., Lang, D., Dalcanton, J. J., et al. 2014, The Astrophysical Journal Supplement Series, 215, 9, doi: 10.1088/0067-0049/215/1/9

[3] Gregersen, D., Seth, A. C., Williams, B. F., et al. 2015, The Astronomical Journal, 150, 189, doi: 10.1088/0004-6256/150/6/189

[4] Williams, B. F., Durbin, M. J., Dalcanton, J. J., et al. 2021, The Astrophysical Journal Supplement Series, 253, 53, doi: 10.3847/1538-4365/abdf4e

[5] Chen, Z., Williams, B., Lang, D., et al. 2024, in American Astronomical Society Meeting Abstracts, Vol. 243, American Astronomical Society Meeting Abstracts, 428.06
