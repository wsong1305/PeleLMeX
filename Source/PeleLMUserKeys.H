/**
  * \brief struct holding PeleLMeX physical BC options
  */
struct BoundaryCondition {
    enum {
        BCInterior = 0,
        BCInflow = 1,
        BCOutflow = 2,
        BCSymmetry = 3,
        BCSlipWallAdiab = 4,
        BCNoSlipWallAdiab = 5,
        BCSlipWallIsotherm = 6,
        BCNoSlipWallIsotherm = 7
    };
};

/**
  * \brief struct holding PeleLMeX physic solver options
           default is LowMachNumber
  */
struct PhysicSolver {
    enum {
        LowMachNumber = 0,
        Incompressible = 1
    };
};

/**
  * \brief struct holding PeleLMeX interpolator options
           default is PiecewiseLinearConserv
  */
struct Interpolator {
    enum {
        PiecewiseLinearConserv = 0,
        PiecewiseLinearConservMinMax = 1,
        PiecewiseConstant = 2
    };
};

/**
  * \brief struct holding PeleLMeX LES SGS model options
           default is None
  */
struct LESModel {
    enum {
        None = 0,
        Smagorinsky = 1,
        WALE = 2,
        Sigma = 3
    };
};

/**
  * \brief struct holding PeleLMeX advection scheme options
           default is Godunov_PLM
  */
struct AdvectionScheme {
    enum {
        Godunov_PLM = 0,
        Godunov_PPM = 1,
        Godunov_PPM_WENOZ = 2,
        Godunov_PPM_NOLIM = 3,
        Godunov_BDS = 4
    };
};

/**
  * \brief struct holding PeleLMeX load balancing cost options
           default is None
  */
struct LoadBalanceCost {
    enum {
        Ncell = 0,
        ChemFunctCallAvg = 1,
        ChemFunctCallMax = 2,
        ChemFunctCallSum = 3,
        UserDefinedDerivedAvg = 4,
        UserDefinedDerivedSum = 5,
        Timers = 6
    };
};

/**
  * \brief struct holding PeleLMeX load balancing method options
           default is SFC
  */
struct LoadBalanceMethod {
    enum {
        SFC = 0,
        Knapsack = 1,
        RoundRobin = 2
    };
};
