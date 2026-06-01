# Single-cell RNA-seq Data

Browse and download the data from the CellxGene portal: https://cellxgene.cziscience.com/e/6f9de485-58cd-4342-bfc4-b3d3dd223aa8.cxg/

# Spatial Transcriptomics Data

We are more than happy to share the spatial transcriptomics data. However, we have not yet identified a suitable platform for easily hosting and sharing these data, particularly the high-resolution H&E images. In contrast, scRNA-seq data have several well-established repositories, such as CELLxGENE and GEO. Therefore, we kindly ask any interested researchers to contact us by email, and we will provide download links upon request. 

Our spatial transcriptomics datasets are organized in the following folder structure. As we are more than happy to share the processed results as much as we can, new contents will be probably added.

## Xenium data

```
xenium
├── img_HE                => H&E images. Raw TIFF files. 
├── img_HE_OME            => H&E images. OME TIFF files. 
├── img_HE_OME_alignment  => Aligning H&E OME images with spatial coordinates. 
├── xenium_data_object    => Data objects in CSV and R RDS files. They include expression matrics, spatial coordinates, cell type/state annotations, niche annotations, etc. 
└── xenium_std_output     => The raw output of the standard Xenium bioinformatics pipeline. 
```

## Visium HD data

```
visiumHD
├── demultiplex_samples     => Dictionary used to split the raw Visium HD data to samples. 
├── visiumhd_data_objects   => Data objects in CSV and R RDS files. They include the expression matrics, spatial cooridnates, high-resolution images, etc. 
├── visiumhd_post_analysis  => Post-analysis results. For example, cell type annotations are included here. 
└── visiumhd_std_output     => The raw output of the standard VisiumHD bioinformatics pipeline. 
```

## Visium

```
visiumST
├── demultiplex_samples    => Dictionary used to split the raw Visium HD data to samples.
├── visiumST_data_objects  => Data objects in CSV and R RDS files. They include the expression matrics, spatial cooridnates, high-resolution images, etc. 
└── visiumST_std_output    => The raw output of the standard Visium bioinformatics pipeline.
```

