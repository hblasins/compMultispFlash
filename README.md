# Computational Multispectral Flash

This repository contains a Matlab implementation of the Computational Multispectral Flash desribed in H. Blasinski and J. Farrell; 'Computational Multispectral Flash,' ICCP 2017.

Please use the following citation if you are using this code for your projects and/or research.
```
@inproceedings{blasinski2017cmf,
    title={Computational Multispectral Flash},
    author={Blasinski, Henryk and Farrell, Joyce},
    booktitle={IEEE International Conference on Computational Photography, ICCP},
    year={2017},
    organization={IEEE}
}
```

## Dependencies

### Required
These additional tools are required if you want to run the code described 
in the paper.

1. [ISET](http://www.imageval.com) - Image Systems Engineering Toolbox. (You can 
request an educational version).
2. [CVX](http://cvxr.com) - An open source Matlab plugin for convex optimization.

### Optional
These tools are necessary if you want to synthesize your own underwater scenes
(or re-synthesize the scenes we used).

3. [RenderToolbox4](http://rendertoolbox.org) - A set of Matlab functions to render images of complex 
scenes. 
4. Coral models - 3D meshes of different corals. We used models from 
[The Hydrous](https://www.thehydro.us), who agreed to share some of their
models for research purposes. We don't own the copyrights, so you will need
to make a separate request.

### Sample data
Sample input data (images) as well as intermediate computation steps can be 
downloaded from the Stanford Digital Repository:
https://purl.stanford.edu/pd601pk8606



