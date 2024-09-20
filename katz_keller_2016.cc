#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <iostream>
#include <string>

// Constants
// set pure component melting points T_m^i at P=0

    std::vector<double> T0;      // Pure component melting points at P=0
    std::vector<double> A;       // Coefficients for linear P-dependence of T_m^i
    std::vector<double> B;       // Coefficients for quadratic P-dependence of T_m^i
    std::vector<double> dS;      // Latent heat of pure components
    std::vector<double> L;       // Latent heat of pure components
    std::string K_T_mode;        // Type of parameterization for K^i(T)
    std::vector<double> r;       // Coefficients for T-dependence of distribution coefficients K^i

    const unsigned int n_components = 3;// stored in some variable


    // already converted to right units
    void
    initialize_parameters() 
    {
    T0 = {1780.0 + 273.15, 1000.0 + 273.15, 710.0 + 273.15, 640.0 + 273.15}; // in K
    A = {45.0e-9, 112.0e-9, 40.8e-9, 30.1e-9};               // K/Pa
    B = {-2.0e-18, -3.37e-18, -1.54e-18, -1.88e-18};         // K/(Pa^2)
    dS = {600e3 / T0[0],
          450e3 / T0[1],
          350e3 / T0[2],
          350e3 / T0[3]};                                    // J/kg/K
    L = {600e3, 450e3, 350e3, 350e3};

    // Initialize r vector (coefficients for T-dependence of distribution coefficients K^i)
    r = {60.0, 30.0, 30.0, 30.0};                            // J/kg/K
    return;
    }


    
    std::vector<double>
    melting_temperatures(const double pressure)
    {
        std::vector<double> Tm (n_components);

        const double Pmax = 6e9;
        if (pressure <= Pmax)
          for (unsigned int i=0; i<n_components; ++i) 
            Tm[i]  =  T0[i] + A[i] * pressure + B[i] * pressure * pressure;

        else
        {
          // safeguard: continue melting point with linear slope above Pmax
          const double dP = 1e7;
          for (unsigned int i=0; i<n_components; ++i) 
          {
            bool ind = pressure > Pmax;
            const double T0_at_Pmax = T0[i] + A[i] * Pmax + B[i] * Pmax * Pmax;
            const double dTdP = ((A[i]*Pmax + B[i]*Pmax*Pmax) - (A[i]*(Pmax-1e7) + B[i]*(Pmax-1e7)*(Pmax-1e7)))/dP;
            Tm[i] = T0_at_Pmax + dTdP * (pressure-Pmax);
          }
        }
        return Tm;
    }



    std::vector<double>
    partition_coefficients(const double pressure, 
                           const double temperature)
    {
        std::vector<double> K (n_components);

        // Implementation in r_DMC is the following instead:
        //std::vector<double> L (n_components);
        //for (unsigned int i=0; i<n_components; ++i) 
        //  L[i] = temperature * dS[i];
        // TODO: ask Tobias Keller which we should use

        std::vector<double> Tm = melting_temperatures(pressure);

        // Parameterization after Rudge, Bercovici, & Spiegelman (2010)
        for (unsigned int i=0; i<n_components; ++i) 
            K[i] = std::exp(L[i]/r[i] * (1./temperature - 1./Tm[i]));

        return K;
    }



    double 
    compute_residual (std::vector<double> composition,
                      std::vector<double> K,
                      bool compute_solidus)
    {
      double residual = -1.;
      for (unsigned int i=0; i<n_components; ++i) 
        if (compute_solidus)
          residual += composition[i]/K[i]; // solidus
        else
          residual += composition[i]*K[i]; // liquidus

      return residual;
    }



//      template <int dim>
      double
//      MeltComponents<dim>::
      T_solidus_liquidus (const double pressure, 
                          std::vector<double> composition, 
                          bool compute_solidus) // const
      {
        // TODO: Exclude invalid compositions (that do not sum up to 1)?

        const std::vector<double> Tm = melting_temperatures(pressure);

        // Set starting guess for Tsol
        const double minTm = *std::min_element(Tm.begin(), Tm.end());
        const double maxTm = *std::max_element(Tm.begin(), Tm.end());
        double mean_Tm = 0.0;

        for (unsigned int i=0; i<n_components; ++i) 
            mean_Tm += composition[i] * Tm[i];

        double T_solidus = std::max(minTm, std::min(maxTm, mean_Tm));

        std::vector<double> K = partition_coefficients(pressure, T_solidus);

        // Get residual for sum(ci_bar/Ki) = 1 or sum(ci_bar*Ki) = 1 (Equations 8 + 9)
        double residual = compute_residual(composition, K, compute_solidus);

        unsigned int n                    =  0;     // initialize iteration count
        const double tolerance            =  1e-10; //1e-10; // tolerance for Newton residual
        const unsigned int max_iterations =  500;   // maximum number of iterations
        const double eps_T                =  5;     // temperature perturbation for finite differencing, degrees

        while (std::abs(residual) > tolerance) 
        {
          // Compute partition coefficients Ki at T+eps_T
          K = partition_coefficients(pressure, T_solidus + eps_T);

          // Get residual at T + eps_T
          double residual_plus_eps_T = compute_residual(composition, K, compute_solidus);

          // Compute partition coefficients Ki at T-eps_T
          K = partition_coefficients(pressure, T_solidus - eps_T);

          // Get residual at T + eps_T
          double residual_minus_eps_T = compute_residual(composition, K, compute_solidus);

          // Finite difference drdT = (r(T+eps_T)-r(T-eps_T))/2/eps_T
          const double dresidualdT  =  (residual_plus_eps_T - residual_minus_eps_T) / (2 * eps_T);

          // Apply Newton correction to current guess of Tsol
          // Note the step size is set to 0.5 whereas the original r_DMC implementation uses 1 
          for (unsigned int i=0; i<n_components; ++i) 
            T_solidus = T_solidus - 0.5 * residual/dresidualdT;

          // Compute partition coefficients Ki at Tsol
          K = partition_coefficients(pressure, T_solidus);

          // Get residual at T_solidus
          residual = compute_residual(composition, K, compute_solidus);

          ++n;

          if (n == max_iterations) 
          {
            std::cerr << "!!! Newton solver for solidus/liquidus T has not converged after " << max_iterations << " iterations !!!" << std::endl;
            break;
          }
        }
      return T_solidus;
    }


int main() 
{
    initialize_parameters();

    const double pressure = 1.5e9; // GPa, as in Fig. B2
    std:: vector<std:: vector<double>> composition = {{1.0, 0.0, 0},  
                                                      {0.95, 0.05, 0},  
                                                      {0.9, 0.1, 0},
                                                      {0.85, 0.15, 0},  
                                                      {0.8, 0.2, 0}, 
                                                      {0.75, 0.25, 0},  
                                                      {0.7, 0.3, 0}, 
                                                      {0.65, 0.35, 0},  
                                                      {0.6, 0.4, 0}, 
                                                      {0.55, 0.45, 0},  
                                                      {0.5, 0.5, 0}, 
                                                      {0.45, 0.55, 0},  
                                                      {0.4, 0.6, 0}, 
                                                      {0.35, 0.65, 0},  
                                                      {0.3, 0.7, 0}, 
                                                      {0.25, 0.75, 0},  
                                                      {0.2, 0.8, 0}, 
                                                      {0.15, 0.85, 0},  
                                                      {0.1, 0.9, 0},  
                                                      {0.05, 0.95, 0}, 
                                                      {0.0, 1.0, 0}}; 

    std::ofstream outputFile("output.txt"); // create a new output file or overwrite an existing one

    if (outputFile.is_open()) 
    { // check if the file was opened successfully

    outputFile << "# Composition Solidus Liquidus" << std::endl;
  
      for (unsigned int j=0; j<composition.size(); ++j)
      {
        outputFile << composition[j][1];
        outputFile << " " << T_solidus_liquidus(pressure, composition[j], true);
        outputFile << " " << T_solidus_liquidus(pressure, composition[j], false) << std::endl;
      }

    outputFile.close(); // close the file when done
    }

    return 0;
}