# Setting Parameters of MCMC dimension reduction

`RiD-kit` uses a JSON-format file to configure MCMC dimension reduction. Here we explain these parameters one by one.

```JSON
{
        "name": "rid-mcmc",
        "init_models": ["model_000.pb","model_001.pb","model_002.pb","model_003.pb"],
        
        "MCMC_Config": {
            "cv_dimension": 2,
            "numb_steps": 10000,
            "numb_walkers": 500,
            "temperature": 300,
            "proj_info": {
                "proj_mode": "cv",
                "proj_cv_index": [0,1]
            },
            "cv_type": "dih",
            "bins": 101
        }
    
},
{
        "name": "rid-mcmc",
        "init_models": ["model_000.pb","model_001.pb","model_002.pb","model_003.pb"],
        
        "MCMC_Config": {
            "cv_dimension": 20,
            "numb_steps": 100000,
            "numb_walkers": 2000,
            "cv_upper_bound": [4.5,5.0,6.0,4.2,4.0,4.0,3.5,4.0,4.0,3.5,3.2,3.0,3.5,3.5,4.0,4.5,4.5,3.75,3.8,4.0],
            "cv_lower_bound": [1.0,2.0,3.5,2.8,2.5,2.0,1.2,1.5,1.5,1.0,1.8,0.5,1.0,1.0,1.5,1.5,1.5,1.75,2.2,2.6],
            "proj_info": {
                "proj_mode": "cv",
                "proj_cv_index": [0,1]
            },
            "cv_type": "dis",
            "bins": 101
        }
    
}
```
* **`"init_models"`** `(List[str])` List of trained free energy models generated in the RiD run.
* **`"cv_dimension"`** `int` An integer to represent the number of CVs used in RiD.
* **`"numb_steps"`** `int` An integer to represent the number of steps for MCMC run.
* **`"numb_walkers"`** `int` Number of parallel walkers used in the MCMC run.
* **`"temperature"`**  `int` Temperature used in the RiD run and MCMC simulation.
* **`"cv_upper_bound"`** `Optional[(List[float])]` Upper bound of the initial value in MCMC for each dimension of CVs, this is usually set for distance CVs.
* **`"cv_lower_bound"`** `Optional[(List[float])]` Lower bound of the initial value in MCMC for each dimension of CVs, this is usually set for distance CVs.
* **`"proj_mode"`** `str` The mode for projecting 2D CV, support `cv` and `path` mode. `cv` mode representing the common projection on the selected CV index, `path` mode representing the projection on the user defined path CV.
* **`"proj_cv_index"`** `(List[int])` Projected index for the 2D CV. Note that 1D projection is done for all CVs.
* **`"cv_type"`** `str` CV type for projection, support `dih` and `dis`, representing dihedral CVs and distance CVs.
* **`"bins"`** `int` Bins for projecting CVs at each dimension.

Another example is given to project 2D path CV.
```JSON
{
        "name": "rid-mcmc",
        "init_models": ["model_000.pb","model_001.pb","model_002.pb","model_003.pb"],
        
        "MCMC_Config": {
            "cv_dimension": 13,
            "numb_steps": 100000,
            "numb_walkers": 2000,
            "cv_upper_bound": [5.33507,4.43402,4.59772,3.40287,4.67268,5.12813,3.28143,3.74967,3.58144,5.34143,7.13014,4.58024,7.41678],
            "cv_lower_bound": [2.38722,0.99005,3.99886,2.08580,3.81799,3.93331,2.768175,3.145618,3.161886,4.661101,5.286805,4.076754,2.893047],
            "proj_info": {
                "proj_mode": "path",
                "proj_cv_index": [0,1,2,3],
                "path_list": [[2.38722,0.99005,3.99886,2.08580],[2.64067,1.20917,3.73831,1.86004],[2.90639,1.50430,3.47935,1.65185],
                [3.18130,1.83919,3.22240,1.46871],[3.46322,2.19574,2.96795,1.32108],[3.75057,2.56491,2.71672,1.22191],
                [4.04218,2.94198,2.46968,1.18344],[4.33720,3.32424,2.22824,1.21147],[4.63498,3.71010,1.99443,1.30170]],
                "path_lm": 0.25
            },
            "cv_type": "dis",
            "bins": 101
        }
    
}
```
In the path CV projection mode, the additonal parameters are:
* **`"proj_cv_index"`** `List[int]` Set of CV index used to defined the path CV.
* **`"path_list"`** `List[float]` List of CVs to defined the path.
* **`"path_lm"`** `float` Lambda in the path CV definition.