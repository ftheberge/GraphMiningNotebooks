import pandas as pd
from bayes_opt import BayesianOptimization
from bayes_opt import UtilityFunction
from h_louvain import hLouvain, load_ABCDH_from_file




class hLouvainBO(hLouvain):

    '''
    Class responsible for Bayesian Optimization made in order to find the best parameters 
    controlling how fast (c) and how much (b) change the coefficient alpha for the hLouvain algorithm

    Parameters
    ----------

    HG : hypernetx.hypergraph.

    hmod_tau : str | float, optional, default= "infinity"
        "infinity" or float greater or equal to 0
        Hyperparameter for hypergraph modularity 
    
    resolution : float, optional, default = 1 (corresponding to the standard settings for modularity calculation) 
        float greater or equal than 0
        Hyperparameter for hypergraph modularity, controlling the weight of degree tax. 
        The classical definition of modularity is retrieved when the resolution parameter is set to 1.

    Note
    ----
    For 'tau', any float >= 0 can be used. 
    Basically it corresponds to hyperedge weights (c/d)^tau when c > d/2 and 0 otherwise.
    Default is "infinity" (strict modularity).
    Other natural choices are 0 -  majority modularity, 1 - linear modularity, 2 - quadratic modularity.

    '''
    
    def set_params(
        self,
        seeds = [1234,5325,5467,4723,999,989,1245, 432,1904,7633],
        xi =1e-3,
        init_points=5,
        n_iter=5,
        pbounds = {'b': (0,1), 'c': (0.01,0.99)},
        bomode = "last_step",
        last_step_top_points = 1,
        show_bo_table = True,
        given_points = []
    ):
        '''
        Function responsible for setting the parameters for hLouvain bayesian optimization procedure
        
        Parameters
        ----------

        seeds : array of integers, optional, default = [1234,5325,5467,4723,999,989,1245, 432,1904,7633]
            The bayesian optimization target function is a mean of hLouvain algoritm executions for seeds taken from this array.
            The length of the array is crucial to control the number of execution.

        xi : float, optional, default = 1e-3
            Theparameter for bayesian optimization, which controls exploration vs exploitation trade off

        init_points : int, optional, default = 5
            The number of initial points (b,c) for bayesian optimization

        n_iter : int, optional, default = 5
            The number of search points for bayesian optimization (bayesian optimization iterations)

        pbounds : dict, default = {'b': (0,1), 'c': (0.01,0.99)}
            Parameter bounds for b and c

        bomode : str, default = "last_step"
            The version of hLouvain bayesian optimization
            The default "last_step" assumes performing the last step optimization only for several best results of core bayesian optimization output.
            Other options: 
            "basic" : no last step, 
            "last_step_all" - performing last step optimization for every execution (not recommended for big networks)
        
        last_step_top_points : int, optional, default = 1
            The parameter used for bomode = "last_step".
            The number of best results for which the last step optimization is performed.
        
        show_bo_table  : boolean, optional, default = True
            When set True the bayesian optimization table is presented

        given_points : list of dicts, optional, default = []
            The points forced to be checked as additional init points for bayesian optimization 
            E.g., given points = [{"b": 0.8, "c": 0.2}, {"b": 0.9, "c": 0.3}]

        '''

        self.seeds = seeds
        self.xi = xi
        self.init_points = init_points
        self.pbounds = pbounds
        self.n_iter = n_iter
        self.bomode = bomode
        self.show_bo_table = show_bo_table
        self.last_step_top_points = last_step_top_points
        #self.given_points = [{"b": 0.8, "c": 0.2}, {"b": 0.9, "c": 0.3}]
        self.given_points = given_points
        self.dts.clear()
        
        self.final_results = pd.DataFrame(
            columns=[
                "b",
                "c",
                "seed",
                "#com",
                "qH",
                "alphas",
                "A",
            ]
        )
     
    
    
    def target_function(self,b,c):
        '''
        The target function for bayesian optimization

        Parameters
        ----------

        b : float
            The parameter from the range (0,1) which controls the how much to change alpha

        c : float
            The parameter from the range (0,1) which controls the how fast to change alpha    
        
        Returns
        -------
        The mean of hypergraph modularity for several execution of hLouvain algorithm with different seeds
        
        '''
        alphas = []
        qHs =[]
        for i in range(30):
            alphas.append(1-((1-b)**i))
        seeds = self.seeds
        for seed in seeds:
            if self.bomode == "basic" or self.bomode == "last_step":
                A, qH, alphas_out = self.h_louvain_community(alphas, 
                                                       change_frequency=c, random_seed=seed)
                if alphas_out[-1] != 1:
                    qH = self.combined_modularity(A, 1, self.hmod_tau, self.resolution)

            if self.bomode == "last_step_all":
                A, A_basic, qH, q_basic, alphas_out = self.h_louvain_community_plus_last_step(alphas = alphas,
                                                       change_frequency=c, random_seed=seed)
                
            
            if len(self.final_results) > 0:
                self.final_results = self.final_results._append(                                    
                                    pd.DataFrame(
                                        [
                                            [b,
                                             c,
                                             seed,
                                             len(A),
                                             qH,
                                             alphas_out,
                                             A,                                                
                                            ]
                                        ],
                                        columns=self.final_results.columns,
                                    ),
                                ignore_index=True,
                            )
            else:
                self.final_results = pd.DataFrame(
                                        [
                                            [b,
                                             c,
                                             seed,
                                             len(A),
                                             qH,
                                             alphas_out,
                                             A,                                                
                                            ]
                                        ],
                                        columns=self.final_results.columns,
                                    )
            
            
            qHs.append(qH)
        result = sum(qHs)/len(seeds)
        return(result)

    def hLouvain_perform_BO(self):

        '''
        The bayesian optimization procedure based on bayesian optimization module


        Returns
        -------

        The dataframe with the final results containing:
        b, c, seed, len(A), qH, alphas_out, A, qh_lstep, A_lstep
        '''    

        optimizer = BayesianOptimization(
            f=self.target_function,
            pbounds=self.pbounds,
            verbose=2, # verbose = 1 prints only when a maximum is observed, verbose = 0 is silent
            random_state=self.random_seed,
        )
        if self.show_bo_table == False:
            optimizer._verbose = 0

        acquisition_function = UtilityFunction(kind="ei", xi=self.xi)
    
        for point in self.given_points:
            optimizer.probe(
                params=point,
                lazy=True,
            )
            
        optimizer.maximize(
            init_points=self.init_points,
            n_iter=self.n_iter,
            acquisition_function=acquisition_function
        )
        
        if self.bomode == "last_step":
            to_check = self.final_results.sort_values(by = ["qH"],  ascending=False).head(self.last_step_top_points)
            qH2_list = []
            A2_list = []
            for index, row in to_check.iterrows():
                A2_item, qH2_item = self.last_step(row["A"], self.hmod_tau,self.resolution)
                A2_list.append(A2_item)
                qH2_list.append(qH2_item)
            to_check = to_check.assign(A_lstep=A2_list)
            to_check = to_check.assign(qH_lstep=qH2_list)
        
            return to_check
        else:
            return self.final_results
            
            



def main():


    HG = load_ABCDH_from_file("datasets/hg-strict_he.txt")


    hBO = hLouvainBO(HG,hmod_tau = "infinity", resolution = 1, random_seed = 875)
    
    hBO.set_params(bomode="last_step", show_bo_table=False)
    hBO_df_results = hBO.hLouvain_perform_BO()

    print(hBO_df_results['qH_lstep'].iloc[0])



if __name__ == "__main__":
    main()
