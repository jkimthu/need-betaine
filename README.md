# need-betaine

## Goal
Image analysis for Deanna's microscopy data of WT vs glycine betaine KO strains.

## Work flow 
1. determine segmentation parameters for different strains (whos_a_cell_dp.m)
2. segment cells from phase images -> gives measurements of "particles" idenfied from images (analysis1.m)
3. plot meaurements with statitics (analysis1.m)

## Initial strategy (analysis1.m)
1. Trim all conditions initially with same segmentation parameters, in effort to avoid human induced bias.
2. Segmentation parameters discriminate between cells and noise.
3. Discrimination is initially genereous. Some non-cell particles may be included, but we will collect far more cells than noise.
	- minWidth = 0.7
	- maxWidth = 3
4. Plot box plots of cell size from each condition
5. After looking at plots, return to image data to confirm conclusions.

### Conclusions from Analysis 1: