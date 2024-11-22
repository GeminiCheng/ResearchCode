# ResearchCode
## Implementation of S3G Module (Social-Support-Based Superpixel Generation)

### Code Description
S3G is the implementation module for the oil spill superpixel segmentation function mentioned in relevant papers.
It is a superpixel generation algorithm based on the SSQM (the Social Support Quantitative Model), using SAM (Spectral Angle Mapper) as the pixel similarity measurement standard and supporting multi parameter input.

The code consists of three parts：
1. Main.m: The main file of S3G, including data reading, superpixel boundaries, result generation, and storage functions.
2. S3G.m: S3G function file, designed to implement superpixel segmentation based on S3G.
3. UEE_evaluation: Evaluation file, used to evaluate the results of S3G based superpixel segmentation. The indicator used is UEE.

### Other information
