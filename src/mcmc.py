import numpy as np
import random

class MCMC(object):
    def __init__(self,n,alpha,beta,jump):
        """
        The value of an MCMC object lies in its simulate method.

        Args:
            n :: Int
            Max number of iterations

            alpha :: Float
            Initial value for step size

            beta :: Float
            Inverse of thermal energy

            jump :: Float -> [Float] -> [Float] -> IO [Float]
            Jumping function (note: must be symmetric about input t)
        Args:
        alpha :: Float
        Step Size
        x :: [?]
        Calibration parameters
        t :: [Float]
        tuning parameters
        Returns:
        t* :: [Float]
        Updated tuning parameters
        """
        self.n 		    = n
        self.alpha     	= alpha
        self.beta       = beta
        self.jump 		= jump

    def simulate(self,x,tInit,simulate,cost):
        """
        Run the MCMC simulation

        Args:

        x :: [?]
        Calibration parameters

        t :: [Float]
        tuning parameters

        simulateFunc :: [Float] -> [Float] -> [Float]
        Runs a simulation. Accepts calibration/tuning parameters.
        Args:
        x :: [?]
        Calibration parameters
        t :: [Float]
        tuning parameters
        Returns:
        y :: [Float]
        Result vector

        costFunc :: [Float] -> [Float] -> [Float] -> Float
        Function for quantifying the accuracy of a TB simulation
        Args:
        x :: [?]
        Calibration parameters
        t :: [Float]
        tuning parameters
        y :: [Float]
        Result vector
        Returns:
            log :: [([Float],Float)]
        A log of the simulation results
        I.e. a length-n list of tuning parameters and associated costs
        """
        # initialize
        # ----------
        log       = [(tInit,float('inf'))] # We WILL update the first time
        best_cost = float('inf')
        for _ in range(self.n):
            t,old_cost = log[-1]
            y 		  = simulate(x,t)
            new_cost  = cost(y)
            deltaCost = new_cost - old_cost
            acc_prob  = np.exp(-self.beta * deltaCost)

            if new_cost < best_cost:
                bestT,best_cost = t,new_cost

            if random.random() < acc_prob:
                log.append( (self.jump(self.alpha, x, t), new_cost) )
            else:
                log.append(log[-1])  # no change

            if random.random() < acc_prob:
                log.append( (self.jump(self.alpha, x, t), new_cost) )
            else:
                log.append(log[-1])  # no change


        return log[1:] #  ignore first point

######################################################################
######################################################################

class SimulateTB(object):
    def __init__(self,pathToReferences,pathToYangTab,pathToTempInput,pathToTempCoords,nEntries,simulateFunc,costFunc):
        """
        A tight-binding simulation.

        Args:
        pathToReferences :: String
        Max number of iterations

        pathToYangTab :: String
        Initial value for step size

        simulateFunc :: [Float] -> [Float] -> [Float]
        Runs TB simulation. Accepts calibration/tuning parameters.
        Args:
        x :: [?]
        Calibration parameters
        t :: [Float]
        tuning parameters
        Returns:
        y :: [Float]
        Result vector

        costFunc :: [Float] -> [Float] -> [Float] -> Float
        Function for quantifying the accuracy of a TB simulation
        Args:
        x :: [?]
        Calibration parameters
        t :: [Float]
        tuning parameters
        y :: [Float]
        Result vector
        Returns:
        cost :: Float
        """

        self.references 		= pathToReferences
        self.yangtab    		= pathToYangTab
        self.pathToTempInput 	= pathToTempInput
        self.pathToTempCoords 	= pathToTempCoords
        self.nEntries 			= nEntries
        self.simulate   		= simulateFunc
        self.cost       		= costFunc
        self.xParams            = (pathToReferences,pathToYangTab,pathToTempInput,pathToTempCoords,nEntries)
    def run_simulation(self,mcmc):
        from parseParam import initializeT
        t0 = initializeT(self.references+'ppots.ref',self.references+'bondints.ref',self.references+'electrons.ref')

        return mcmc.simulate(self.xParams,t0,self.simulate,self.cost)


######################################################################
# Functions as arguments to MCMC/TBparam objects
######################################################################

def defaultCostFunc(y):
    eDFT,eDFT_sd,eTB,fDFT,fDFT_sd,fTB = y

    eDiff       = eDFT - eTB
    eScalar     = np.dot(eDiff,eDiff)/len(eTB)
    forceDiff   = fDFT-fTB
    forceScalar = np.dot(forceDiff,forceDiff)/len(fTB)

    return eScalar/eDFT_sd + forceScalar/fDFT_sd


def CosineSimilarityCostFunc(y):
    eDFT,eDFT_sd,eTB,fDFT,fDFT_sd,fTB = y
    similarity_E = np.dot(eDFT,eTB) /(np.linalg.norm(eDFT) * np.linalg.norm(eTB))
    similarity_F = np.dot(fDFT,fTB) /(np.linalg.norm(fDFT) * np.linalg.norm(fTB))
    return (1 - similarity_E) + (1 - similarity_F)

def varyAllJumpFunc(alpha,x,t):
    sdev = 1 * alpha 						  # replace the 1 with a function of x,t? Do we really want all parameters to vary uniformly?
    changeVec = np.multiply(t,np.random.normal(scale=sdev)) # centering at input t ensures symmetric distribution
    return t + changeVec

def varyOneJumpFunc(alpha,x,t):
    sdev = 1 * alpha 						  # replace the 1 with a function of x,t? Do we really want all parameters to vary uniformly?
    changeInd = random.choice(range(len(t)))
    t[changeInd]+= sdev * np.random.normal(1) * t[changeInd]
    return t



######################################################################
# Other functions
######################################################################
def analyzeSensitivity(log):
    # Given a list of t's (and costs, if that's useful)
    # need to evaluate sensitivity of each tuning parameter
    # simple solution: standard deviation of parameters?

    log = log[300:] # remove 'burn-in' samples?
    for t,c in log: raise NotImplementedError


######################################################################
######################################################################
def test():
    """
    Simulate the TB model using default parameters for TB model and MCMC

    THINGS WE WANT TO TEST
    MCMC params
    Vary nIter(1000-100000),alpha(0.01 - 1),beta(1-1.001)
    TB params
    Just stick with hydrocarbons at first, then move on to O and N
    Cost Function
    just consider forces? just consider energies?
    least squared differences normalized by standard deviation?
    cosine similarity cost function?
    Jumping distribution function
    Vary one parameter at a time vs all parameters at once?
    Analysis of results
    Various ways of quantifying how much each param changes (TBD)

    """
    from inputLATTE import tightBind

    latteRoot = '/LATTE/'
    TB_1   = SimulateTB(latteRoot+'TBparam',latteRoot+'yangtab_24'
                            ,latteRoot+'bl/inputblock.dat',latteRoot+'TBparam'
                            ,10,tightBind,defaultCostFunc)
    MCMC_1 = MCMC(10,1,10,varyAllJumpFunc)

    output =  TB_1.run_simulation(MCMC_1)

    print output
    # analyze(output)
