====================================
README for '3DS soil model' (soil3ds)
====================================

This is '3DS soil model' (soil3ds), a 3D Soil model of soil for water and N balances adapted from the STICS soil module.

See 
- Louarn G, Faverjon L, Migault V, Escobar-Gutiérrez A, Combes D. (2016). Assessment of ‘3DS’, a soil module for individual-based models of plant communities. In: IEEE International Conference on Functional-Structural Plant Growth Modeling, Simulation, Visualization and Applications (FSPMA), 125–132. doi: 10.1109/FSPMA.2016.7818298
- N. Brisson, M. Launay, B. Mary, N. Beaudoin (2008) Conceptual basis, formalisations and parameterization of the STICS crop model, Quae, Versailles.
- N. Brisson & A. Perrier (1991). A semi-empirical model of bare soil evaporation for crop simulation models. Water resources research, 27(5), 719-727. 
- Lebon, E., Dumas, V., Pieri, P., & Schultz, H. R. (2003). Modelling the seasonal dynamics of the soil water balance of vineyards. Functional Plant Biology, 30(6), 699-710.


## 1. Getting Started

These instructions will get you a copy of *soil3ds* up and running on your local 
machine.

### 1.1 Prerequisites

To install and use *soil3ds*, you need first to install the dependencies.

*soil3ds* has been tested on Windows.
 
#### 1.1.1 Install the dependencies on Windows 10 64 bit
1) Create a conda environment with miniconda3
    ```bash
    conda create -n myenvname python=3.7 xlrd=2.0.1 numpy=1.20.3 scipy=1.7.3 pandas=1.3.4
    ```

2) Place yourself in the created environment  : `conda activate myenvname`

3) Install *soil3ds*
    1) Git console :
        ```bash
        git clone https://github.com/glouarn/soil3ds.git
        ```
    2) installation in the conda environment (in folder `soil3ds`)
        ```bash
        python setup.py develop
        ```


### 1.3 Running

To run a simulation example :

* 1. place yourself in folder `soil3ds/test`
  2. run from the console:
		```bash
        python test_soil_coupling1.py
        ```


## 2. Reading the docs

To build the user and reference guides:


## 3. Testing

The test allows to verify that the model implementation accurately 
represents the developer’s conceptual description of the model and its solution.


## Contact

For any question, send an email to <gaetan.louarn@inrae.fr>.

## Authors

**Gaëtan LOUARN**, **Eric LEBON** - see file [AUTHORS](AUTHORS) for details

## License

This project is licensed under the CeCILL-C License - see file [LICENSE](LICENSE) for details
